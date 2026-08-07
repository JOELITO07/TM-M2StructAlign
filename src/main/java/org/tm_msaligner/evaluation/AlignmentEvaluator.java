package org.tm_msaligner.evaluation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;

/**
 * Computes reference-based residue-pair metrics and a topology-aware gap rate.
 * Residues are identified by protein ID and ungapped sequence position, so the
 * two MSAs may have different lengths and gap patterns.
 */
public final class AlignmentEvaluator {
  private final FastaAlignment preparedReference;
  private final TopologyDataset preparedTopology;
  private final Path preparedReferenceFile;
  private final List<String> preparedProteinIds;
  private final PairSets preparedReferencePairs;

  /** Creates a stateless evaluator for one-off comparisons. */
  public AlignmentEvaluator() {
    preparedReference = null;
    preparedTopology = null;
    preparedReferenceFile = null;
    preparedProteinIds = null;
    preparedReferencePairs = null;
  }

  /**
   * Prepares a reusable evaluator. The reference alignment and topology are
   * parsed and indexed only once, which is preferable for seed directories.
   */
  public AlignmentEvaluator(Path referenceFile, Path topologyFile) throws IOException {
    preparedReference = FastaAlignment.read(referenceFile);
    preparedTopology = TopologyDataset.read(topologyFile);
    preparedReferenceFile = referenceFile;
    preparedProteinIds = new ArrayList<>(preparedReference.sequences().keySet());
    validateReference(preparedReference, preparedTopology, referenceFile, topologyFile);
    preparedReferencePairs = collectPairs(
        preparedReference, preparedTopology, preparedProteinIds);
  }

  public AlignmentEvaluationResult evaluate(
      Path referenceFile,
      Path testFile,
      Path topologyFile) throws IOException {

    return new AlignmentEvaluator(referenceFile, topologyFile).evaluate(testFile);
  }

  /** Evaluates one MSA using the reference supplied to the constructor. */
  public AlignmentEvaluationResult evaluate(Path testFile) throws IOException {
    if (preparedReference == null) {
      throw new IllegalStateException(
          "This evaluator has no prepared reference. Use AlignmentEvaluator(reference, topology) "
              + "or evaluate(reference, test, topology).");
    }

    FastaAlignment test = FastaAlignment.read(testFile);

    validateTest(
        preparedReference,
        test,
        preparedReferenceFile,
        testFile);

    PairSets referencePairs = preparedReferencePairs;
    PairSets testPairs = collectPairs(test, preparedTopology, preparedProteinIds);

    long correctPairs = intersectionSize(referencePairs.allPairs(), testPairs.allPairs());
    long correctTmPairs = intersectionSize(referencePairs.tmPairs(), testPairs.tmPairs());

    double pairPrecision = ratio(
        correctPairs,
        testPairs.allPairs().size(),
        referencePairs.allPairs().isEmpty());
    double pairRecall = ratio(
        correctPairs,
        referencePairs.allPairs().size(),
        testPairs.allPairs().isEmpty());

    double tmPairPrecision = ratio(
        correctTmPairs,
        testPairs.tmPairs().size(),
        referencePairs.tmPairs().isEmpty());
    double tmPairRecall = ratio(
        correctTmPairs,
        referencePairs.tmPairs().size(),
        testPairs.tmPairs().isEmpty());

    GapCounts tmGaps = countTmGaps(test, preparedTopology, preparedProteinIds);

    return new AlignmentEvaluationResult(
        preparedProteinIds.size(),
        referencePairs.allPairs().size(),
        testPairs.allPairs().size(),
        correctPairs,
        pairPrecision,
        pairRecall,
        harmonicMean(pairPrecision, pairRecall),
        referencePairs.tmPairs().size(),
        testPairs.tmPairs().size(),
        correctTmPairs,
        tmPairPrecision,
        tmPairRecall,
        harmonicMean(tmPairPrecision, tmPairRecall),
        tmGaps.gaps(),
        tmGaps.opportunities(),
        tmGaps.opportunities() == 0
            ? Double.NaN
            : (double) tmGaps.gaps() / tmGaps.opportunities());
  }

  private void validateReference(
      FastaAlignment reference,
      TopologyDataset topology,
      Path referenceFile,
      Path topologyFile) throws IOException {
    for (String id : reference.sequences().keySet()) {
      AAArray annotatedResidues = topology.residues(id);
      if (annotatedResidues == null) {
        throw new IOException(
            "Protein '" + id + "' is missing from topology file " + topologyFile);
      }

      String expected = residueSequence(annotatedResidues);
      String referenceUngapped = removeGaps(reference.sequences().get(id));
      if (!expected.equals(referenceUngapped)) {
        throw new IOException(
            "Ungapped sequence mismatch for protein '" + id + "' between topology "
                + topologyFile + " and reference " + referenceFile
                + ". topologyLength=" + expected.length()
                + ", referenceLength=" + referenceUngapped.length());
      }
    }
  }

  private void validateTest(
      FastaAlignment reference,
      FastaAlignment test,
      Path referenceFile,
      Path testFile) throws IOException {

    Set<String> referenceIds = reference.sequences().keySet();
    Set<String> testIds = test.sequences().keySet();
    if (!referenceIds.equals(testIds)) {
      Set<String> missing = new HashSet<>(referenceIds);
      missing.removeAll(testIds);
      Set<String> unexpected = new HashSet<>(testIds);
      unexpected.removeAll(referenceIds);
      throw new IOException(
          "Protein IDs differ between reference " + referenceFile + " and test " + testFile
              + ". Missing=" + missing + ", unexpected=" + unexpected);
    }

    for (String id : referenceIds) {
      String referenceUngapped = removeGaps(reference.sequences().get(id));
      String testUngapped = removeGaps(test.sequences().get(id));
      if (!referenceUngapped.equals(testUngapped)) {
        throw new IOException(
            "Ungapped sequence mismatch for protein '" + id + "' between reference "
                + referenceFile + " and test " + testFile
                + ". referenceLength=" + referenceUngapped.length()
                + ", testLength=" + testUngapped.length());
      }
    }
  }

  private PairSets collectPairs(
      FastaAlignment alignment,
      TopologyDataset topology,
      List<String> proteinIds) {

    Set<ResiduePair> allPairs = new HashSet<>();
    Set<ResiduePair> tmPairs = new HashSet<>();
    int[] residuePositions = new int[proteinIds.size()];

    for (int column = 0; column < alignment.alignmentLength(); column++) {
      List<AlignedResidue> residues = new ArrayList<>();
      for (int sequenceIndex = 0; sequenceIndex < proteinIds.size(); sequenceIndex++) {
        String id = proteinIds.get(sequenceIndex);
        char letter = alignment.sequences().get(id).charAt(column);
        if (!isGap(letter)) {
          int residueIndex = residuePositions[sequenceIndex]++;
          AA residue = topology.residues(id).AAAt(residueIndex);
          residues.add(new AlignedResidue(id, residueIndex, residue.getType().isTMRegion()));
        }
      }

      for (int first = 0; first < residues.size() - 1; first++) {
        for (int second = first + 1; second < residues.size(); second++) {
          AlignedResidue left = residues.get(first);
          AlignedResidue right = residues.get(second);
          ResiduePair pair = new ResiduePair(
              left.proteinId(), left.position(), right.proteinId(), right.position());
          allPairs.add(pair);
          if (left.transmembrane() && right.transmembrane()) {
            tmPairs.add(pair);
          }
        }
      }
    }
    return new PairSets(allPairs, tmPairs);
  }

  private GapCounts countTmGaps(
      FastaAlignment alignment,
      TopologyDataset topology,
      List<String> proteinIds) {

    int[] residuePositions = new int[proteinIds.size()];
    long tmGapCount = 0;
    long tmResidueCount = 0;

    for (String id : proteinIds) {
      AAArray residues = topology.residues(id);
      for (int residueIndex = 0; residueIndex < residues.getSize(); residueIndex++) {
        AA residue = residues.AAAt(residueIndex);
        if (residue.getType().isTMRegion()) {
          tmResidueCount++;
        }
      }
    }

    for (int column = 0; column < alignment.alignmentLength(); column++) {
      boolean[] gaps = new boolean[proteinIds.size()];
      boolean[] transmembrane = new boolean[proteinIds.size()];

      for (int sequenceIndex = 0; sequenceIndex < proteinIds.size(); sequenceIndex++) {
        String id = proteinIds.get(sequenceIndex);
        char letter = alignment.sequences().get(id).charAt(column);
        gaps[sequenceIndex] = isGap(letter);
        if (!gaps[sequenceIndex]) {
          AA residue = topology.residues(id).AAAt(residuePositions[sequenceIndex]++);
          transmembrane[sequenceIndex] = residue.getType().isTMRegion();
        }
      }

      for (int first = 0; first < proteinIds.size() - 1; first++) {
        for (int second = first + 1; second < proteinIds.size(); second++) {
          if (gaps[first] && !gaps[second] && transmembrane[second]) {
            tmGapCount++;
          } else if (!gaps[first] && transmembrane[first] && gaps[second]) {
            tmGapCount++;
          }
        }
      }
    }

    long opportunities = tmResidueCount * Math.max(0, proteinIds.size() - 1L);
    return new GapCounts(tmGapCount, opportunities);
  }

  private static long intersectionSize(Set<ResiduePair> first, Set<ResiduePair> second) {
    Set<ResiduePair> smaller = first.size() <= second.size() ? first : second;
    Set<ResiduePair> larger = first.size() <= second.size() ? second : first;
    long count = 0;
    for (ResiduePair pair : smaller) {
      if (larger.contains(pair)) {
        count++;
      }
    }
    return count;
  }

  private static double ratio(long numerator, long denominator, boolean otherSetIsEmpty) {
    if (denominator == 0) {
      return otherSetIsEmpty ? 1.0 : 0.0;
    }
    return (double) numerator / denominator;
  }

  private static double harmonicMean(double precision, double recall) {
    return precision + recall == 0.0
        ? 0.0
        : 2.0 * precision * recall / (precision + recall);
  }

  private static String removeGaps(String sequence) {
    return sequence.replace("-", "").replace(".", "").toUpperCase();
  }

  private static String residueSequence(AAArray residues) {
    StringBuilder sequence = new StringBuilder(residues.getSize());
    for (int index = 0; index < residues.getSize(); index++) {
      sequence.append(Character.toUpperCase(residues.AAAt(index).getLetter()));
    }
    return sequence.toString();
  }

  private static boolean isGap(char letter) {
    return letter == AA.GAP_IDENTIFIER || letter == '.';
  }

  private record AlignedResidue(String proteinId, int position, boolean transmembrane) {
  }

  private record ResiduePair(
      String firstProtein,
      int firstPosition,
      String secondProtein,
      int secondPosition) {
  }

  private record PairSets(Set<ResiduePair> allPairs, Set<ResiduePair> tmPairs) {
  }

  private record GapCounts(long gaps, long opportunities) {
  }
}
