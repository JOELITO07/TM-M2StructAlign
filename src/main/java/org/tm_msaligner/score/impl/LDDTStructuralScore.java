package org.tm_msaligner.score.impl;

import java.util.Arrays;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.structure.StructuralNeighborhood;
import org.tm_msaligner.util.AA;

/**
 * Symmetric, coverage-aware and sparse structural consistency score inspired by
 * FoldMason MSTA lDDT.
 *
 * <p>Local neighbours inside {@code R0} are precomputed once from each dense
 * distance matrix and reused by every solution evaluation. Because the current
 * dataset contains only {@code *_D.csv} and {@code *_idx.csv}, the contact
 * weights are derived from the distance matrix instead of being loaded from a
 * separate {@code *_W.csv} file:</p>
 *
 * <pre>
 * weight(i,j) = 1 - D(i,j) / R0, for D(i,j) &lt; R0
 * </pre>
 *
 * <p>The same query-contact weight contributes to the numerator and the
 * denominator. Missing, gapped or structurally unmapped target neighbours
 * remain in the denominator and contribute zero to the numerator.</p>
 */
public class LDDTStructuralScore implements StructuralScore {

  private static final float MIN_CONTACT_WEIGHT = 1.0e-6f;

  private final float R0;
  private volatile StructuralCache structuralCache;

  public LDDTStructuralScore(float R0) {
    if (!Float.isFinite(R0) || R0 <= 0.0f) {
      throw new IllegalArgumentException("R0 must be finite and greater than zero");
    }
    this.R0 = R0;
  }

  @Override
  public <S extends StructuralTM_MSASolution> double compute(S solution, AA[][] msa) {
    int numberOfSequences = solution.variables().size();
    int alignmentLength = solution.getAlignmentLength();

    if (numberOfSequences < 2 || alignmentLength == 0) {
      return 0.0;
    }

    StructuralCache cache = getOrBuildStructuralCache(solution, numberOfSequences);

    double[] columnScoreSum = new double[alignmentLength];
    int[] columnScoreCount = new int[alignmentLength];
    int[] residuesPerColumn = new int[alignmentLength];
    int[][] structToColumn = new int[numberOfSequences][];

    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      StructuralNeighborhood neighborhood = cache.neighborhoods[sequence];
      structToColumn[sequence] = neighborhood == null
          ? new int[0]
          : buildStructToColumnMap(msa[sequence], neighborhood.residueCount());

      for (int column = 0; column < alignmentLength; column++) {
        if (!msa[sequence][column].isGap()) {
          residuesPerColumn[column]++;
        }
      }
    }

    // Ordered pairs evaluate A->B and B->A, preserving sequence-order symmetry.
    for (int query = 0; query < numberOfSequences; query++) {
      StructuralNeighborhood queryNeighborhood = cache.neighborhoods[query];
      if (queryNeighborhood == null) {
        continue;
      }

      for (int target = 0; target < numberOfSequences; target++) {
        if (query == target || cache.distanceMatrices[target] == null) {
          continue;
        }

        float[][] targetDistances = cache.distanceMatrices[target];

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
              queryNeighborhood,
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

    // Columns without structural support contribute zero.
    return total / alignmentLength;
  }

  /**
   * Builds all sparse neighbourhoods in two passes and caches them for reuse.
   * The cache is rebuilt only if a different set of distance-matrix objects is
   * supplied.
   */
  private synchronized StructuralCache getOrBuildStructuralCache(
      StructuralTM_MSASolution solution,
      int numberOfSequences) {

    float[][][] currentMatrices = new float[numberOfSequences][][];
    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      currentMatrices[sequence] = solution.getDistanceMatrixByIndex(sequence);
    }

    StructuralCache currentCache = structuralCache;
    if (currentCache != null && currentCache.matches(currentMatrices)) {
      return currentCache;
    }

    StructuralNeighborhood[] neighborhoods = new StructuralNeighborhood[numberOfSequences];
    for (int sequence = 0; sequence < numberOfSequences; sequence++) {
      float[][] distances = currentMatrices[sequence];
      neighborhoods[sequence] = distances == null ? null : buildNeighborhood(distances);
    }

    StructuralCache rebuilt = new StructuralCache(currentMatrices, neighborhoods);
    structuralCache = rebuilt;
    return rebuilt;
  }

  /**
   * Constructs a sparse neighbourhood from one dense distance matrix.
   *
   * <p>First pass: count valid neighbours inside R0 for every residue.</p>
   * <p>Second pass: allocate exact-size arrays and populate indices, distances
   * and distance-derived weights.</p>
   */
  private StructuralNeighborhood buildNeighborhood(float[][] distances) {
    int residueCount = distances.length;
    int[] neighborCounts = new int[residueCount];

    // First pass: count neighbours.
    for (int residue = 0; residue < residueCount; residue++) {
      if (distances[residue] == null || distances[residue].length != residueCount) {
        throw new IllegalArgumentException("Distance matrix must be square");
      }

      for (int neighbor = 0; neighbor < residueCount; neighbor++) {
        if (neighbor == residue) {
          continue;
        }

        float distance = distances[residue][neighbor];
        if (Float.isFinite(distance) && distance >= 0.0f && distance < R0) {
          neighborCounts[residue]++;
        }
      }
    }

    int[][] neighborIndices = new int[residueCount][];
    float[][] neighborDistances = new float[residueCount][];
    float[][] neighborWeights = new float[residueCount][];

    for (int residue = 0; residue < residueCount; residue++) {
      int count = neighborCounts[residue];
      neighborIndices[residue] = new int[count];
      neighborDistances[residue] = new float[count];
      neighborWeights[residue] = new float[count];
    }

    // Second pass: populate exact-size arrays.
    int[] insertionPositions = new int[residueCount];
    for (int residue = 0; residue < residueCount; residue++) {
      for (int neighbor = 0; neighbor < residueCount; neighbor++) {
        if (neighbor == residue) {
          continue;
        }

        float distance = distances[residue][neighbor];
        if (!Float.isFinite(distance) || distance < 0.0f || distance >= R0) {
          continue;
        }

        int position = insertionPositions[residue]++;
        neighborIndices[residue][position] = neighbor;
        neighborDistances[residue][position] = distance;
        neighborWeights[residue][position] = contactWeight(distance);
      }
    }

    return new StructuralNeighborhood(
        residueCount,
        neighborIndices,
        neighborDistances,
        neighborWeights);
  }

  /**
   * Computes weighted query-to-target lDDT for one aligned anchor pair.
   * The sparse query neighbourhood defines the denominator.
   */
  private double directionalLddt(
      StructuralNeighborhood queryNeighborhood,
      float[][] targetDistances,
      AA[] targetAlignment,
      int[] queryStructToColumn,
      int queryAnchorIndex,
      int targetAnchorIndex) {

    if (queryAnchorIndex >= queryNeighborhood.residueCount()
        || targetAnchorIndex >= targetDistances.length
        || targetDistances[targetAnchorIndex] == null) {
      return Double.NaN;
    }

    int[] neighborIndices = queryNeighborhood.neighborIndices(queryAnchorIndex);
    float[] neighborDistances = queryNeighborhood.neighborDistances(queryAnchorIndex);
    float[] neighborWeights = queryNeighborhood.neighborWeights(queryAnchorIndex);

    if (neighborIndices.length == 0) {
      return Double.NaN;
    }

    double weightedScoreSum = 0.0;
    double weightSum = 0.0;

    for (int position = 0; position < neighborIndices.length; position++) {
      int queryNeighborIndex = neighborIndices[position];
      float queryWeight = neighborWeights[position];

      // Every valid query contact remains in the denominator.
      weightSum += queryWeight;

      if (queryNeighborIndex >= queryStructToColumn.length) {
        continue;
      }

      int neighborColumn = queryStructToColumn[queryNeighborIndex];
      if (neighborColumn < 0 || neighborColumn >= targetAlignment.length) {
        continue;
      }

      AA targetNeighbor = targetAlignment[neighborColumn];
      if (targetNeighbor.isGap() || targetNeighbor.getStructIndex() < 0) {
        continue;
      }

      int targetNeighborIndex = targetNeighbor.getStructIndex();
      if (targetNeighborIndex >= targetDistances[targetAnchorIndex].length) {
        continue;
      }

      float targetDistance = targetDistances[targetAnchorIndex][targetNeighborIndex];
      if (!Float.isFinite(targetDistance)) {
        continue;
      }

      float difference = Math.abs(neighborDistances[position] - targetDistance);
      weightedScoreSum += queryWeight * thresholdScore(difference);
    }

    return weightSum > 0.0 ? weightedScoreSum / weightSum : Double.NaN;
  }

  /**
   * Distance-decay contact weight derived directly from the available D matrix.
   * No external W matrix is required.
   */
  private float contactWeight(float distance) {
    return Math.max(MIN_CONTACT_WEIGHT, 1.0f - distance / R0);
  }

  /**
   * Equivalent to testing the four lDDT thresholds 0.5, 1, 2 and 4 Angstrom,
   * without allocating or iterating over a threshold array.
   */
  private double thresholdScore(float difference) {
    if (difference >= 4.0f) {
      return 0.0;
    } else if (difference >= 2.0f) {
      return 0.25;
    } else if (difference >= 1.0f) {
      return 0.50;
    } else if (difference >= 0.5f) {
      return 0.75;
    }
    return 1.0;
  }

  private int[] buildStructToColumnMap(AA[] alignedSequence, int residueCount) {
    int[] map = new int[residueCount];
    Arrays.fill(map, -1);

    for (int column = 0; column < alignedSequence.length; column++) {
      AA aa = alignedSequence[column];
      int structIndex = aa.isGap() ? -1 : aa.getStructIndex();
      if (structIndex >= 0 && structIndex < residueCount && map[structIndex] < 0) {
        map[structIndex] = column;
      }
    }

    return map;
  }

  private static final class StructuralCache {
    private final float[][][] distanceMatrices;
    private final StructuralNeighborhood[] neighborhoods;

    private StructuralCache(
        float[][][] distanceMatrices,
        StructuralNeighborhood[] neighborhoods) {
      this.distanceMatrices = distanceMatrices;
      this.neighborhoods = neighborhoods;
    }

    private boolean matches(float[][][] otherMatrices) {
      if (distanceMatrices.length != otherMatrices.length) {
        return false;
      }
      for (int sequence = 0; sequence < distanceMatrices.length; sequence++) {
        if (distanceMatrices[sequence] != otherMatrices[sequence]) {
          return false;
        }
      }
      return true;
    }
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
    return "Symmetric sparse weighted and coverage-aware reference-free MSTA lDDT";
  }

  public String getDescription() {
    return description();
  }

  @Override
  public String getName() {
    return name();
  }
}
