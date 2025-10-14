package org.tm_msaligner.util.structure;

import java.util.Collections;
import java.util.Map;

/**
 * Container for structural descriptors extracted from an AlphaFold2 model. The structure stores a
 * residue-wise mapping that can be queried from different parts of the algorithm.
 */
public class AlphaFoldStructure {
  private final Map<Integer, ResidueStructuralInfo> residueInfo;
  private final double[][] distanceMatrix;
  private final double contactThreshold;

  public AlphaFoldStructure(
      Map<Integer, ResidueStructuralInfo> residueInfo,
      double[][] distanceMatrix,
      double contactThreshold) {
    this.residueInfo = Map.copyOf(residueInfo);
    this.distanceMatrix = distanceMatrix;
    this.contactThreshold = contactThreshold;
  }

  public Map<Integer, ResidueStructuralInfo> getResidueInfo() {
    return Collections.unmodifiableMap(residueInfo);
  }

  public ResidueStructuralInfo getResidue(int residueIndex) {
    return residueInfo.get(residueIndex);
  }

  public double[][] getDistanceMatrix() {
    return distanceMatrix;
  }

  public double getContactThreshold() {
    return contactThreshold;
  }
}
