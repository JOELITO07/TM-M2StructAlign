package org.tm_msaligner.util.structure;

/**
 * Enumeration representing the coarse-grained secondary structure classes that are extracted from
 * AlphaFold2 models. The type is intentionally simple to make it easy to connect the structural
 * descriptors with the sequence alignment objective functions.
 */
public enum SecondaryStructureType {
  HELIX,
  STRAND,
  COIL;

  /**
   * Maps a character representing the secondary structure (as found in PDB files) to the
   * corresponding enumeration value.
   *
   * <p>By default, characters associated with helices (H, G, I) are mapped to {@link #HELIX},
   * characters associated with extended beta conformations (E, B) are mapped to {@link #STRAND},
   * and all other characters map to {@link #COIL}.
   *
   * @param code secondary structure code. If {@code null} the method returns {@link #COIL}.
   * @return the corresponding {@link SecondaryStructureType}.
   */
  public static SecondaryStructureType fromCode(Character code) {
    if (code == null) {
      return COIL;
    }

    return switch (Character.toUpperCase(code)) {
      case 'H', 'G', 'I' -> HELIX;
      case 'E', 'B' -> STRAND;
      default -> COIL;
    };
  }
}
