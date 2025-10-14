package org.tm_msaligner.score.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import org.tm_msaligner.score.Score;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.structure.AlphaFoldStructure;
import org.tm_msaligner.util.structure.ResidueContact;
import org.tm_msaligner.util.structure.ResidueStructuralInfo;
import org.tm_msaligner.util.structure.SecondaryStructureType;

/**
 * Score that favours alignments respecting AlphaFold2-derived structural constraints. The metric is
 * composed of two parts: (i) a contact preservation score that rewards alignments where residues in
 * contact remain contiguous (i.e. without gaps inserted between them) and (ii) a secondary
 * structure agreement component that penalises placing structured residues (helices and strands)
 * against incompatible structural contexts.
 */
public class StructurePreservingScore implements Score {

  private static final double DEFAULT_CONTACT_WEIGHT = 0.6;
  private static final double DEFAULT_SECONDARY_STRUCTURE_WEIGHT = 0.4;

  private final Map<String, AlphaFoldStructure> structuresBySequence;
  private final double contactWeight;
  private final double secondaryStructureWeight;

  /**
   * Creates the score with default weights for the contact and secondary structure components.
   *
   * @param structures mapping between sequence identifiers and their AlphaFold2-derived
   *     structural descriptors. Keys are normalised internally (trimmed and converted to upper
   *     case) to ease integration with FASTA headers used across the project.
   */
  public StructurePreservingScore(Map<String, AlphaFoldStructure> structures) {
    this(structures, DEFAULT_CONTACT_WEIGHT, DEFAULT_SECONDARY_STRUCTURE_WEIGHT);
  }

  /**
   * Creates the score using custom weights for the contact and secondary structure terms.
   *
   * @param structures mapping between sequence identifiers and AlphaFold2 structures. Keys are
   *     normalised internally (trimmed and converted to upper case).
   * @param contactWeight weight assigned to the contact-preserving term. Must be non-negative.
   * @param secondaryStructureWeight weight assigned to the secondary-structure agreement term.
   *     Must be non-negative. At least one of the weights must be positive.
   */
  public StructurePreservingScore(
      Map<String, AlphaFoldStructure> structures,
      double contactWeight,
      double secondaryStructureWeight) {
    Objects.requireNonNull(structures, "structures");
    if (structures.isEmpty()) {
      throw new IllegalArgumentException("At least one AlphaFold structure must be provided");
    }
    if (contactWeight < 0 || secondaryStructureWeight < 0) {
      throw new IllegalArgumentException("Weights must be non-negative");
    }
    if (contactWeight + secondaryStructureWeight == 0) {
      throw new IllegalArgumentException("At least one weight must be positive");
    }

    this.structuresBySequence = normaliseKeys(structures);
    this.contactWeight = contactWeight;
    this.secondaryStructureWeight = secondaryStructureWeight;
  }

  @Override
  public String getName() {
    return "StructurePreservingScore";
  }

  @Override
  public String getDescription() {
    return "Multi-objective score that rewards AlphaFold2-guided structural coherence";
  }

  @Override
  public boolean isAMinimizationScore() {
    return false;
  }

  @Override
  public String name() {
    return getName();
  }

  @Override
  public String description() {
    return getDescription();
  }

  @Override
  public <S extends TM_MSASolution> double compute(S solution, AA[][] decodedSequences) {
    Objects.requireNonNull(solution, "solution");
    Objects.requireNonNull(decodedSequences, "decodedSequences");
    if (decodedSequences.length == 0) {
      return 0.0;
    }

    List<SequenceStructuralContext> contexts = buildContexts(solution, decodedSequences);
    if (contexts.isEmpty()) {
      // No structural information is available for the sequences involved in the solution. In this
      // case the score returns a neutral value so that it does not bias the optimisation process.
      return 1.0;
    }

    double contactScore = computeContactPreservationScore(contexts);
    double secondaryScore = computeSecondaryStructureScore(contexts);

    double totalWeight = contactWeight + secondaryStructureWeight;
    double weightedContact = contactWeight * contactScore;
    double weightedSecondary = secondaryStructureWeight * secondaryScore;

    return (weightedContact + weightedSecondary) / totalWeight;
  }

  private List<SequenceStructuralContext> buildContexts(
      TM_MSASolution solution, AA[][] decodedSequences) {
    List<SequenceStructuralContext> contexts = new ArrayList<>();

    for (int seqIndex = 0; seqIndex < decodedSequences.length; seqIndex++) {
      AlphaFoldStructure structure = resolveStructureForSequence(solution, seqIndex);
      if (structure == null) {
        continue;
      }

      SecondaryStructureType[] columnTypes = new SecondaryStructureType[decodedSequences[seqIndex].length];
      Map<Integer, Integer> residueToColumn =
          mapResiduesToColumns(decodedSequences[seqIndex], structure, columnTypes);
      contexts.add(
          new SequenceStructuralContext(
              structure, residueToColumn, columnTypes, decodedSequences[seqIndex]));
    }

    return contexts;
  }

  private AlphaFoldStructure resolveStructureForSequence(TM_MSASolution solution, int seqIndex) {
    Object attribute = solution.getAttribute("SeqName" + seqIndex);
    if (attribute == null) {
      return null;
    }

    String name = attribute.toString();
    return findStructureForName(name);
  }

  private AlphaFoldStructure findStructureForName(String name) {
    String normalised = normaliseKey(name);
    AlphaFoldStructure structure = structuresBySequence.get(normalised);
    if (structure != null) {
      return structure;
    }

    int pipeIndex = normalised.indexOf('|');
    if (pipeIndex >= 0) {
      structure = structuresBySequence.get(normalised.substring(pipeIndex + 1).trim());
      if (structure != null) {
        return structure;
      }
      structure = structuresBySequence.get(normalised.substring(0, pipeIndex).trim());
      if (structure != null) {
        return structure;
      }
    }

    int spaceIndex = normalised.indexOf(' ');
    if (spaceIndex >= 0) {
      String truncated = normalised.substring(0, spaceIndex).trim();
      structure = structuresBySequence.get(truncated);
      if (structure != null) {
        return structure;
      }
    }

    return null;
  }

  private Map<Integer, Integer> mapResiduesToColumns(
      AA[] alignedSequence, AlphaFoldStructure structure, SecondaryStructureType[] columnTypes) {
    Map<Integer, Integer> mapping = new HashMap<>();
    List<Integer> sortedResidues = new ArrayList<>(structure.getResidueInfo().keySet());
    sortedResidues.sort(Comparator.naturalOrder());

    int residueCursor = 0;
    for (int column = 0; column < alignedSequence.length && residueCursor < sortedResidues.size(); column++) {
      if (alignedSequence[column].isGap()) {
        continue;
      }

      int residueIndex = sortedResidues.get(residueCursor++);
      mapping.put(residueIndex, column);
      ResidueStructuralInfo info = structure.getResidue(residueIndex);
      columnTypes[column] = info != null ? info.getSecondaryStructure() : null;
    }

    return mapping;
  }

  private double computeContactPreservationScore(List<SequenceStructuralContext> contexts) {
    double weightedScoreSum = 0.0;
    double weightAccumulator = 0.0;

    for (SequenceStructuralContext context : contexts) {
      AlphaFoldStructure structure = context.structure();
      double sequenceScoreSum = 0.0;
      double sequenceWeightSum = 0.0;

      for (ResidueStructuralInfo info : structure.getResidueInfo().values()) {
        Integer columnI = context.residueToColumn().get(info.getResidueIndex());
        if (columnI == null) {
          continue;
        }
        for (ResidueContact contact : info.getContacts()) {
          if (info.getResidueIndex() >= contact.getResidueIndex()) {
            continue; // Count each pair once.
          }

          Integer columnJ = context.residueToColumn().get(contact.getResidueIndex());
          if (columnJ == null) {
            continue;
          }

          double distance = contact.getDistance();
          double pairWeight = Double.isNaN(distance)
              ? 1.0
              : Math.max(0.0, 1.0 - (distance / structure.getContactThreshold()));
          if (pairWeight == 0.0) {
            pairWeight = 1.0;
          }

          double gapPenalty = computeGapPenalty(context.sequence(), columnI, columnJ);
          sequenceScoreSum += gapPenalty * pairWeight;
          sequenceWeightSum += pairWeight;
        }
      }

      if (sequenceWeightSum > 0) {
        weightedScoreSum += sequenceScoreSum / sequenceWeightSum;
        weightAccumulator += 1.0;
      }
    }

    return weightAccumulator > 0 ? weightedScoreSum / weightAccumulator : 0.0;
  }

  private double computeGapPenalty(AA[] sequence, int columnA, int columnB) {
    int start = Math.min(columnA, columnB);
    int end = Math.max(columnA, columnB);
    int gapCount = 0;

    for (int column = start; column <= end; column++) {
      if (sequence[column].isGap()) {
        gapCount++;
      }
    }

    return 1.0 / (1.0 + gapCount);
  }

  private double computeSecondaryStructureScore(List<SequenceStructuralContext> contexts) {
    if (contexts.size() < 2) {
      return 1.0;
    }

    int alignmentLength = contexts.get(0).columnTypes().length;
    double scoreSum = 0.0;
    int evaluatedColumns = 0;

    for (int column = 0; column < alignmentLength; column++) {
      double columnScore = evaluateColumnSecondaryStructure(contexts, column);
      if (!Double.isNaN(columnScore)) {
        scoreSum += columnScore;
        evaluatedColumns++;
      }
    }

    return evaluatedColumns > 0 ? scoreSum / evaluatedColumns : 0.0;
  }

  private double evaluateColumnSecondaryStructure(
      List<SequenceStructuralContext> contexts, int column) {
    double scoreSum = 0.0;
    int comparisons = 0;

    for (int i = 0; i < contexts.size(); i++) {
      SecondaryStructureType typeI = contexts.get(i).columnTypes()[column];
      if (typeI == null || typeI == SecondaryStructureType.COIL) {
        continue;
      }
      for (int j = i + 1; j < contexts.size(); j++) {
        SecondaryStructureType typeJ = contexts.get(j).columnTypes()[column];
        if (typeJ == null || typeJ == SecondaryStructureType.COIL) {
          continue;
        }

        comparisons++;
        scoreSum += (typeI == typeJ) ? 1.0 : 0.0;
      }
    }

    if (comparisons == 0) {
      return Double.NaN;
    }
    return scoreSum / comparisons;
  }

  private Map<String, AlphaFoldStructure> normaliseKeys(Map<String, AlphaFoldStructure> input) {
    Map<String, AlphaFoldStructure> normalised = new HashMap<>();
    for (Map.Entry<String, AlphaFoldStructure> entry : input.entrySet()) {
      if (entry.getValue() == null) {
        continue;
      }
      String key = normaliseKey(entry.getKey());
      normalised.put(key, entry.getValue());
    }
    return Collections.unmodifiableMap(normalised);
  }

  private String normaliseKey(String key) {
    return Optional.ofNullable(key).map(value -> value.trim().toUpperCase()).orElse("");
  }

  private record SequenceStructuralContext(
      AlphaFoldStructure structure,
      Map<Integer, Integer> residueToColumn,
      SecondaryStructureType[] columnTypes,
      AA[] sequence) {}
}
