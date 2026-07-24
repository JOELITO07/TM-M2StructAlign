package org.tm_msaligner.score.impl;

import java.util.HashMap;
import java.util.Map;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;

/**
 * Reference-free structural consistency score inspired by FoldMason MSTA lDDT.
 *
 * <p>The score is evaluated in both directions for every structure pair. Every
 * structural neighbour of the query residue inside {@code R0} remains in the
 * denominator, even when that neighbour is not aligned in the target. Pairwise
 * residue scores are mapped back to MSA columns, combined with column occupancy
 * through a harmonic mean, and normalized by the complete alignment length.</p>
 */
public class LDDTStructuralScore implements StructuralScore {

  private final float R0;
  private final float[] thresholds = new float[] {0.5f, 1.0f, 2.0f, 4.0f};

  public LDDTStructuralScore(float R0) {
    this.R0 = R0;
  }

  @Override
  public <S extends StructuralTM_MSASolution> double compute(S solution, AA[][] msa) {
    int numberOfSequences = solution.variables().size();
    int alignmentLength = solution.getAlignmentLength();

    if (numberOfSequences < 2 || alignmentLength == 0) {
      return 0.0;
    }

    double[] columnScoreSum = new double[alignmentLength];
    int[] columnScoreCount = new int[alignmentLength];
    int[] residuesPerColumn = new int[alignmentLength];


    Map<Integer, Integer>[] structToColumn = new Map[numberOfSequences];

    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      structToColumn[sequence] = buildStructToColumnMap(msa[sequence]);
      for (int column = 0; column < alignmentLength; column++) {
        if (!msa[sequence][column].isGap()) {
          residuesPerColumn[column]++;
        }
      }
    }

    // Ordered pairs evaluate A->B and B->A, making the final score symmetric.
    for (int query = 0; query < numberOfSequences; query++) {
      float[][] queryDistances = solution.getDistanceMatrixByIndex(query);
      if (queryDistances == null) {
        continue;
      }

      for (int target = 0; target < numberOfSequences; target++) {
        if (query == target) {
          continue;
        }

        float[][] targetDistances = solution.getDistanceMatrixByIndex(target);
        if (targetDistances == null) {
          continue;
        }

        for (int column = 0; column < alignmentLength; column++) {
          AA queryAnchor = msa[query][column];
          AA targetAnchor = msa[target][column];

          if (queryAnchor.isGap() || targetAnchor.isGap()) {
            continue;
          }

          int queryAnchorIndex = queryAnchor.getStructIndex();
          int targetAnchorIndex = targetAnchor.getStructIndex();
          if (queryAnchorIndex < 0 || targetAnchorIndex < 0) {
            continue;
          }

          double score = directionalLddt(
              queryDistances,
              targetDistances,
              msa[target],
              structToColumn[query],
              queryAnchorIndex,
              targetAnchorIndex);

          if (!Double.isNaN(score)) {
            columnScoreSum[column] += score;
            columnScoreCount[column]++;
          }
        }
      }
    }

    double total = 0.0;
    for (int column = 0; column < alignmentLength; column++) {
      double occupancy = residuesPerColumn[column] / (double) numberOfSequences;
      double structuralScore = columnScoreCount[column] == 0
          ? 0.0
          : columnScoreSum[column] / columnScoreCount[column];

      if (occupancy > 0.0 && structuralScore > 0.0) {
        total += 2.0 * structuralScore * occupancy / (structuralScore + occupancy);
      }
    }

    // Columns without structural support contribute zero, preserving coverage sensitivity.
    return total / alignmentLength;
  }

  /**
   * Computes query-to-target lDDT for one aligned anchor pair.
   * All query neighbours within R0 are included in the denominator. A missing,
   * gapped, or structurally unmapped target neighbour contributes zero.
   */
  private double directionalLddt(
      float[][] queryDistances,
      float[][] targetDistances,
      AA[] targetAlignment,
      Map<Integer, Integer> queryStructToColumn,
      int queryAnchorIndex,
      int targetAnchorIndex) {

    if (queryAnchorIndex >= queryDistances.length
        || targetAnchorIndex >= targetDistances.length) {
      return Double.NaN;
    }

    double scoreSum = 0.0;
    int neighbourCount = 0;

    for (int queryNeighbourIndex = 0;
         queryNeighbourIndex < queryDistances[queryAnchorIndex].length;
         queryNeighbourIndex++) {

      if (queryNeighbourIndex == queryAnchorIndex) {
        continue;
      }

      float queryDistance = queryDistances[queryAnchorIndex][queryNeighbourIndex];
      if (!Float.isFinite(queryDistance) || queryDistance >= R0) {
        continue;
      }

      // The complete query neighbourhood defines the denominator.
      neighbourCount++;

      Integer neighbourColumn = queryStructToColumn.get(queryNeighbourIndex);
      if (neighbourColumn == null
          || neighbourColumn < 0
          || neighbourColumn >= targetAlignment.length) {
        continue;
      }

      AA targetNeighbour = targetAlignment[neighbourColumn];
      if (targetNeighbour.isGap() || targetNeighbour.getStructIndex() < 0) {
        continue;
      }

      int targetNeighbourIndex = targetNeighbour.getStructIndex();
      if (targetNeighbourIndex >= targetDistances[targetAnchorIndex].length) {
        continue;
      }

      float targetDistance = targetDistances[targetAnchorIndex][targetNeighbourIndex];
      if (!Float.isFinite(targetDistance)) {
        continue;
      }

      scoreSum += thresholdScore(Math.abs(queryDistance - targetDistance));
    }

    return neighbourCount == 0 ? Double.NaN : scoreSum / neighbourCount;
  }

  private double thresholdScore(double difference) {
    int preserved = 0;
    for (float threshold : thresholds) {
      if (difference < threshold) {
        preserved++;
      }
    }
    return preserved / (double) thresholds.length;
  }

  private Map<Integer, Integer> buildStructToColumnMap(AA[] alignedSequence) {
    Map<Integer, Integer> map = new HashMap<>();
    for (int column = 0; column < alignedSequence.length; column++) {
      AA aa = alignedSequence[column];
      if (!aa.isGap() && aa.getStructIndex() >= 0) {
        map.putIfAbsent(aa.getStructIndex(), column);
      }
    }
    return map;
  }

  @Override
  public boolean isAMinimizationScore() {
    return false;
  }

  @Override
  public String name() {
    return "lDDT_structural";
  }

  @Override
  public String description() {
    return "Symmetric, coverage-aware reference-free MSTA lDDT structural consistency";
  }

  public String getDescription() {
    return description();
  }

  @Override
  public String getName() {
    return name();
  }
}
