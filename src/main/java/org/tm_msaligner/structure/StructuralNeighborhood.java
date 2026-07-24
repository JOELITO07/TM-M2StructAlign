package org.tm_msaligner.structure;

/**
 * Immutable sparse representation of the local structural neighbourhood of a
 * protein. For every structural residue, only neighbours inside the configured
 * distance cutoff are stored.
 */
public final class StructuralNeighborhood {

  private final int residueCount;
  private final int[][] neighborIndices;
  private final float[][] neighborDistances;
  private final float[][] neighborWeights;

  public StructuralNeighborhood(
      int residueCount,
      int[][] neighborIndices,
      float[][] neighborDistances,
      float[][] neighborWeights) {

    if (residueCount < 0) {
      throw new IllegalArgumentException("residueCount must be non-negative");
    }
    if (neighborIndices == null
        || neighborDistances == null
        || neighborWeights == null
        || neighborIndices.length != residueCount
        || neighborDistances.length != residueCount
        || neighborWeights.length != residueCount) {
      throw new IllegalArgumentException("Invalid structural-neighborhood dimensions");
    }

    for (int residue = 0; residue < residueCount; residue++) {
      int length = neighborIndices[residue].length;
      if (neighborDistances[residue].length != length
          || neighborWeights[residue].length != length) {
        throw new IllegalArgumentException(
            "Neighbor indices, distances and weights must have equal lengths");
      }
    }

    this.residueCount = residueCount;
    this.neighborIndices = neighborIndices;
    this.neighborDistances = neighborDistances;
    this.neighborWeights = neighborWeights;
  }

  public int residueCount() {
    return residueCount;
  }

  public int[] neighborIndices(int residueIndex) {
    return neighborIndices[residueIndex];
  }

  public float[] neighborDistances(int residueIndex) {
    return neighborDistances[residueIndex];
  }

  public float[] neighborWeights(int residueIndex) {
    return neighborWeights[residueIndex];
  }

  public int neighborCount(int residueIndex) {
    return neighborIndices[residueIndex].length;
  }
}
