package org.tm_msaligner.util.structure;

/**
 * Represents a contact between two residues derived from an AlphaFold2 model. Contacts are
 * computed using the distance between the C-alpha atoms of the residues.
 */
public class ResidueContact {
  private final int residueIndex;
  private final double distance;

  public ResidueContact(int residueIndex, double distance) {
    this.residueIndex = residueIndex;
    this.distance = distance;
  }

  public int getResidueIndex() {
    return residueIndex;
  }

  public double getDistance() {
    return distance;
  }
}
