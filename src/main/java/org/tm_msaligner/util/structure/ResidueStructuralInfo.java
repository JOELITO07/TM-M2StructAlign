package org.tm_msaligner.util.structure;

import java.util.Collections;
import java.util.List;

/**
 * Stores residue-level descriptors derived from an AlphaFold2 structure. The information can be
 * used as constraints or additional objectives when evaluating sequence alignments.
 */
public class ResidueStructuralInfo {
  private final int residueIndex;
  private final String residueName;
  private final SecondaryStructureType secondaryStructure;
  private final double confidence;
  private final List<ResidueContact> contacts;

  public ResidueStructuralInfo(
      int residueIndex,
      String residueName,
      SecondaryStructureType secondaryStructure,
      double confidence,
      List<ResidueContact> contacts) {
    this.residueIndex = residueIndex;
    this.residueName = residueName;
    this.secondaryStructure = secondaryStructure;
    this.confidence = confidence;
    this.contacts = List.copyOf(contacts);
  }

  public int getResidueIndex() {
    return residueIndex;
  }

  public String getResidueName() {
    return residueName;
  }

  public SecondaryStructureType getSecondaryStructure() {
    return secondaryStructure;
  }

  public double getConfidence() {
    return confidence;
  }

  public List<ResidueContact> getContacts() {
    return Collections.unmodifiableList(contacts);
  }
}
