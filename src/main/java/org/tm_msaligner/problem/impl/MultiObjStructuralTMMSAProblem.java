package org.tm_msaligner.problem.impl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.tm_msaligner.problem.StructuralTMMSAProblem;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;
import org.uma.jmetal.util.errorchecking.JMetalException;

public class MultiObjStructuralTMMSAProblem extends StructuralTMMSAProblem {

  private final List<StructuralScore> scoreList;

  public MultiObjStructuralTMMSAProblem(
      String msaProblemFileName,
      List<StructuralScore> scoreList,
      List<String> preComputedFiles,
      String distanceDatasetDir,
      String name) throws IOException {
    super(msaProblemFileName, preComputedFiles, distanceDatasetDir);

    setNumberOfObjectives(scoreList.size());
    setName(name);
    this.scoreList = scoreList;
  }

  /**
   * Reads a precomputed alignment and restores the canonical input order by
   * matching normalized FASTA identifiers. A trailing .pdb extension is
   * removed from both input and precomputed identifiers, e.g.
   * ada2b_human.pdb -> ada2b_human.
   */
  @Override
  public List<AAArray> readDataFromFastaFile(String dataFile) throws IOException {
    LinkedHashMap<String, ProteinSequence> fastaSequences =
        FastaReaderHelper.readFastaProteinSequence(new File(dataFile));

    Map<String, AAArray> sequencesById = new HashMap<>();
    for (Map.Entry<String, ProteinSequence> entry : fastaSequences.entrySet()) {
      ProteinSequence protein = entry.getValue();
      String header = protein.getOriginalHeader();
      if (header == null || header.isBlank()) {
        header = entry.getKey();
      }

      String normalizedId = normalizeProteinId(header);
      AAArray previous = sequencesById.put(
          normalizedId,
          new AAArray(protein.getSequenceAsString()));

      if (previous != null) {
        throw new JMetalException(
            "Duplicated normalized protein ID '" + normalizedId + "' in " + dataFile);
      }
    }

    List<AAArray> orderedSequences = new ArrayList<>(listOfSequenceNames.size());
    int expectedAlignmentLength = -1;

    for (int index = 0; index < listOfSequenceNames.size(); index++) {
      String expectedId = normalizeProteinId(listOfSequenceNames.get(index).toString());
      AAArray alignedSequence = sequencesById.remove(expectedId);

      if (alignedSequence == null) {
        throw new JMetalException(
            "Protein '" + expectedId + "' was not found in precomputed alignment " + dataFile);
      }

      validateUngappedSequence(
          expectedId,
          originalSequences.get(index),
          alignedSequence,
          dataFile);

      if (expectedAlignmentLength < 0) {
        expectedAlignmentLength = alignedSequence.getSize();
      } else if (alignedSequence.getSize() != expectedAlignmentLength) {
        throw new JMetalException(
            "Different aligned sequence lengths in " + dataFile
                + ": expected " + expectedAlignmentLength
                + " but protein '" + expectedId + "' has " + alignedSequence.getSize());
      }

      orderedSequences.add(alignedSequence);
    }

    if (!sequencesById.isEmpty()) {
      throw new JMetalException(
          "Unexpected proteins in " + dataFile + ": " + sequencesById.keySet());
    }

    return orderedSequences;
  }

  private String normalizeProteinId(String header) {
    if (header == null) {
      throw new JMetalException("Null FASTA header");
    }

    String value = header.trim();
    if (value.startsWith(">")) {
      value = value.substring(1).trim();
    }

    if (value.isEmpty()) {
      throw new JMetalException("Empty FASTA identifier");
    }

    // Keep only the identifier token before an optional description.
    value = value.split("\\s+", 2)[0];

    String[] fields = value.split("\\|");
    if (fields.length >= 2
        && ("sp".equalsIgnoreCase(fields[0]) || "tr".equalsIgnoreCase(fields[0]))) {
      value = fields[1];
    } else if (fields.length >= 1) {
      value = fields[0];
    }

    value = value.trim();
    if (value.toLowerCase().endsWith(".pdb")) {
      value = value.substring(0, value.length() - 4);
    }

    if (value.isEmpty()) {
      throw new JMetalException("Empty normalized FASTA identifier derived from header: " + header);
    }

    return value.toLowerCase();
  }

  private void validateUngappedSequence(
      String proteinId,
      AAArray expected,
      AAArray aligned,
      String dataFile) {

    String expectedSequence = expected.toString().replace("-", "").toUpperCase();
    String alignedSequence = aligned.toString().replace("-", "").toUpperCase();

    if (!expectedSequence.equals(alignedSequence)) {
      throw new JMetalException(
          "Ungapped sequence mismatch for protein '" + proteinId + "' in " + dataFile
              + ". Expected length=" + expectedSequence.length()
              + ", obtained length=" + alignedSequence.length());
    }
  }

  @Override
  public StructuralTM_MSASolution evaluate(StructuralTM_MSASolution solution) {
    solution.removeGapColumns();
    AA[][] decodedSequences = solution.decodeToMatrix();

    for (int i = 0; i < numberOfObjectives(); i++) {
      solution.objectives()[i] = scoreList.get(i).compute(solution, decodedSequences)
          * (scoreList.get(i).isAMinimizationScore() ? 1.0 : -1.0);
    }

    return solution;
  }
}
