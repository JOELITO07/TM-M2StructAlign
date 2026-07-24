package org.tm_msaligner.crossover;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.util.errorchecking.Check;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/**
 * Topology-safe single-point crossover for structural TM multiple alignments.
 *
 * <p>The cut is defined as a residue boundary, not merely as an aligned column.
 * A candidate is accepted only when the equivalent boundary is valid in both
 * parents and does not split a continuous TM helix in any sequence. When no
 * safe boundary exists, the corresponding parent is copied unchanged.</p>
 */
public class BioSPXMSACrossover implements CrossoverOperator<StructuralTM_MSASolution> {
  private final JMetalRandom randomGenerator;
  private final double probability;

  /*
   * Kept for constructor/API compatibility with existing runners. Unsafe cuts
   * are no longer accepted probabilistically, regardless of alpha.
   */
  private final double alpha;

  public BioSPXMSACrossover(double probability, double alpha) {
    Check.probabilityIsValid(probability);
    Check.probabilityIsValid(alpha);

    this.randomGenerator = JMetalRandom.getInstance();
    this.probability = probability;
    this.alpha = alpha;
  }

  @Override
  public List<StructuralTM_MSASolution> execute(List<StructuralTM_MSASolution> parents) {
    Check.notNull(parents);
    Check.that(parents.size() == 2, "The number of parents is not 2");
    Check.that(
        parents.get(0).variables().size() == parents.get(1).variables().size(),
        "The two parents have different number of sequences: "
            + parents.get(0).variables().size()
            + ", "
            + parents.get(1).variables().size());

    return doCrossover(parents);
  }

  private List<StructuralTM_MSASolution> doCrossover(
      List<StructuralTM_MSASolution> parents) {

    StructuralTM_MSASolution parent1 = parents.get(0);
    StructuralTM_MSASolution parent2 = parents.get(1);

    List<StructuralTM_MSASolution> children = new ArrayList<>(2);
    children.add(msaCrossover(parent1, parent2));
    children.add(msaCrossover(parent2, parent1));

    return children;
  }

  private StructuralTM_MSASolution msaCrossover(
      StructuralTM_MSASolution parentA,
      StructuralTM_MSASolution parentB) {

    if (randomGenerator.nextDouble() >= probability) {
      return new StructuralTM_MSASolution(parentA);
    }

    CutPlan cutPlan = selectCompatibleCut(parentA, parentB);
    if (cutPlan == null) {
      return new StructuralTM_MSASolution(parentA);
    }

    int cut = cutPlan.cutInParentA;
    List<List<Integer>> gapsGroupFirstBlock = new ArrayList<>();
    List<Integer> gapsGroup;

    /*
     * Copy the prefix gap groups from parentA up to and including the selected
     * aligned column.
     */
    for (int sequence = 0; sequence < parentA.variables().size(); sequence++) {
      gapsGroup = parentA.variables().get(sequence);

      List<Integer> gaps = new ArrayList<>();
      int numberOfPrefixGaps = 0;

      for (int group = 1; group < gapsGroup.size(); group += 2) {
        if (cut >= gapsGroup.get(group)) {
          gaps.add(gapsGroup.get(group - 1));
          gaps.add(gapsGroup.get(group));
          numberOfPrefixGaps += gapsGroup.get(group) - gapsGroup.get(group - 1) + 1;
        } else {
          if (cut >= gapsGroup.get(group - 1)) {
            gaps.add(gapsGroup.get(group - 1));
            gaps.add(cut);
            numberOfPrefixGaps += cut - gapsGroup.get(group - 1) + 1;
          }
          break;
        }
      }

      int calculatedResiduesOnLeft = cut - numberOfPrefixGaps + 1;
      if (calculatedResiduesOnLeft != cutPlan.residuesOnLeft[sequence]) {
        throw new IllegalStateException(
            "Inconsistent crossover residue mapping for sequence "
                + sequence
                + ": calculated="
                + calculatedResiduesOnLeft
                + ", planned="
                + cutPlan.residuesOnLeft[sequence]);
      }

      gapsGroupFirstBlock.add(gaps);
    }

    List<Integer> positionsToCutParentB = new ArrayList<>(
        cutPlan.suffixStartInParentB.length);
    for (int position : cutPlan.suffixStartInParentB) {
      positionsToCutParentB.add(position);
    }

    int minimumSuffixPosition = Collections.min(positionsToCutParentB);

    /*
     * Append the suffix gap groups from parentB. The padding step keeps all
     * sequence suffixes vertically synchronized when their equivalent boundary
     * occurs at different aligned columns.
     */
    for (int sequence = 0; sequence < parentB.variables().size(); sequence++) {
      int suffixPosition = positionsToCutParentB.get(sequence);
      int outputPosition = cut;
      List<Integer> gaps = gapsGroupFirstBlock.get(sequence);

      if (suffixPosition > minimumSuffixPosition) {
        int shift = suffixPosition - minimumSuffixPosition;

        if (!gaps.isEmpty()) {
          int lastGapEnd = gaps.get(gaps.size() - 1);
          if (lastGapEnd != outputPosition) {
            gaps.add(outputPosition + 1);
            gaps.add(outputPosition + shift);
          } else {
            gaps.set(gaps.size() - 1, outputPosition + shift);
          }
        } else {
          gaps.add(outputPosition + 1);
          gaps.add(outputPosition + shift);
        }

        outputPosition += shift;
      }

      gapsGroup = parentB.variables().get(sequence);
      for (int group = 0; group < gapsGroup.size(); group += 2) {
        if (gapsGroup.get(group) >= suffixPosition) {
          gaps.add(outputPosition + (gapsGroup.get(group) - suffixPosition) + 1);
          gaps.add(outputPosition + (gapsGroup.get(group + 1) - suffixPosition) + 1);
        }
      }
    }

    StructuralTM_MSASolution child = new StructuralTM_MSASolution(
        parentA.getStructuralTMMSAProblem(),
        gapsGroupFirstBlock);

    child.mergeGapsGroups();
    return child;
  }

  /**
   * Returns every aligned column in parentA that represents a valid residue
   * boundary in both parents and does not split a TM helix.
   *
   * <p>Package-private for focused unit testing.</p>
   */
  List<Integer> findCompatibleCutColumns(
      StructuralTM_MSASolution parentA,
      StructuralTM_MSASolution parentB) {

    int numberOfSequences = parentA.variables().size();
    if (numberOfSequences != parentB.variables().size()) {
      return Collections.emptyList();
    }

    int alignmentLength = parentA.getAlignmentLength();
    if (alignmentLength < 3) {
      return Collections.emptyList();
    }

    int[][] prefixResidueCounts = buildPrefixResidueCounts(parentA);
    int[][] suffixStartLookup = buildSuffixStartLookup(parentB);
    List<Integer> compatibleCuts = new ArrayList<>();

    /*
     * The cut is the last column supplied by parentA. Exclude terminal
     * boundaries so that both prefix and suffix contain aligned columns.
     */
    for (int cut = 1; cut < alignmentLength - 1; cut++) {
      if (parentA.isGapColumn(cut)) {
        continue;
      }

      if (createCutPlan(
          parentA,
          parentB,
          cut,
          prefixResidueCounts,
          suffixStartLookup) != null) {
        compatibleCuts.add(cut);
      }
    }

    return compatibleCuts;
  }

  /**
   * Tests one explicit cut column. Package-private for unit tests.
   */
  boolean isCompatibleCut(
      StructuralTM_MSASolution parentA,
      StructuralTM_MSASolution parentB,
      int cut) {

    if (parentA.variables().size() != parentB.variables().size()
        || cut < 1
        || cut >= parentA.getAlignmentLength() - 1
        || parentA.isGapColumn(cut)) {
      return false;
    }

    return createCutPlan(
        parentA,
        parentB,
        cut,
        buildPrefixResidueCounts(parentA),
        buildSuffixStartLookup(parentB)) != null;
  }

  private CutPlan selectCompatibleCut(
      StructuralTM_MSASolution parentA,
      StructuralTM_MSASolution parentB) {

    List<Integer> compatibleCuts = findCompatibleCutColumns(parentA, parentB);
    if (compatibleCuts.isEmpty()) {
      return null;
    }

    int selectedCut = compatibleCuts.get(
        randomGenerator.nextInt(0, compatibleCuts.size() - 1));

    return createCutPlan(
        parentA,
        parentB,
        selectedCut,
        buildPrefixResidueCounts(parentA),
        buildSuffixStartLookup(parentB));
  }

  private CutPlan createCutPlan(
      StructuralTM_MSASolution parentA,
      StructuralTM_MSASolution parentB,
      int cut,
      int[][] prefixResidueCounts,
      int[][] suffixStartLookup) {

    int numberOfSequences = parentA.variables().size();
    int[] residuesOnLeft = new int[numberOfSequences];
    int[] suffixStartInParentB = new int[numberOfSequences];

    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      int leftResidueCount = prefixResidueCounts[sequence][cut];
      residuesOnLeft[sequence] = leftResidueCount;

      /*
       * Validate the biological boundary in both parents. Both parents normally
       * share the same original sequences, but checking both makes the operator
       * contract explicit and protects against future data-source changes.
       */
      if (splitsTMHelixAfterResidue(parentA, sequence, leftResidueCount)
          || splitsTMHelixAfterResidue(parentB, sequence, leftResidueCount)) {
        return null;
      }

      if (leftResidueCount < 0
          || leftResidueCount >= suffixStartLookup[sequence].length) {
        return null;
      }

      int suffixStart = suffixStartLookup[sequence][leftResidueCount];
      if (suffixStart < 0
          || suffixStart > parentB.getAlignmentLength(sequence)) {
        return null;
      }

      suffixStartInParentB[sequence] = suffixStart;
    }

    return new CutPlan(cut, residuesOnLeft, suffixStartInParentB);
  }

  /**
   * Cumulative number of non-gap residues at every aligned column.
   */
  private int[][] buildPrefixResidueCounts(StructuralTM_MSASolution solution) {
    int numberOfSequences = solution.variables().size();
    int alignmentLength = solution.getAlignmentLength();
    int[][] prefixCounts = new int[numberOfSequences][alignmentLength];

    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      int residueCount = 0;
      for (int column = 0; column < alignmentLength; column++) {
        if (solution.getSeqPos1BasedAtColumn(sequence, column) > 0) {
          residueCount++;
        }
        prefixCounts[sequence][column] = residueCount;
      }
    }

    return prefixCounts;
  }

  /**
   * Maps a residue count r to the aligned column where parentB's suffix starts
   * after retaining exactly r residues in the prefix.
   */
  private int[][] buildSuffixStartLookup(StructuralTM_MSASolution solution) {
    int numberOfSequences = solution.variables().size();
    int[][] suffixStart = new int[numberOfSequences][];

    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      int originalLength = solution.getOriginalSequences().get(sequence).getSize();
      suffixStart[sequence] = new int[originalLength + 1];
      Arrays.fill(suffixStart[sequence], -1);
      suffixStart[sequence][0] = 0;

      int alignmentLength = solution.getAlignmentLength(sequence);
      for (int column = 0; column < alignmentLength; column++) {
        int residuePosition = solution.getSeqPos1BasedAtColumn(sequence, column);
        if (residuePosition > 0 && residuePosition <= originalLength) {
          suffixStart[sequence][residuePosition] = column + 1;
        }
      }
    }

    return suffixStart;
  }

  /**
   * Returns true when the boundary after leftResidueCount residues lies between
   * two consecutive TM residues in the original ungapped sequence.
   *
   * <p>Package-private for focused unit testing.</p>
   */
  boolean splitsTMHelixAfterResidue(
      StructuralTM_MSASolution solution,
      int sequenceIndex,
      int leftResidueCount) {

    AAArray original = solution.getOriginalSequences().get(sequenceIndex);
    if (leftResidueCount <= 0 || leftResidueCount >= original.getSize()) {
      return false;
    }

    AA leftResidue = original.AAAt(leftResidueCount - 1);
    AA rightResidue = original.AAAt(leftResidueCount);

    return isTM(leftResidue) && isTM(rightResidue);
  }

  private boolean isTM(AA residue) {
    return residue != null
        && residue.getType() != null
        && residue.getType().isTMRegion();
  }

  private static final class CutPlan {
    private final int cutInParentA;
    private final int[] residuesOnLeft;
    private final int[] suffixStartInParentB;

    private CutPlan(
        int cutInParentA,
        int[] residuesOnLeft,
        int[] suffixStartInParentB) {

      this.cutInParentA = cutInParentA;
      this.residuesOnLeft = residuesOnLeft;
      this.suffixStartInParentB = suffixStartInParentB;
    }
  }

  @Override
  public int numberOfRequiredParents() {
    return 2;
  }

  @Override
  public int numberOfGeneratedChildren() {
    return 2;
  }

  @Override
  public double crossoverProbability() {
    return probability;
  }

  /**
   * Returns the retained compatibility parameter. It no longer relaxes the
   * strict TM-boundary constraint.
   */
  public double alpha() {
    return alpha;
  }
}
