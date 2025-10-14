package org.tm_msaligner.util.structure;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;
import org.biojava.nbio.structure.secstruc.SecStrucType;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.GroupType;

/**
 * Utility capable of parsing structural information from AlphaFold2 PDB files. The parser extracts
 * per-residue descriptors (secondary structure, confidence scores, and contacts) that can later be
 * used to bias or constrain the multiple sequence alignment process.
 */
public class AlphaFoldStructureParser {
  private static final double DEFAULT_CONTACT_THRESHOLD = 8.0;

  private final double contactThreshold;

  public AlphaFoldStructureParser() {
    this(DEFAULT_CONTACT_THRESHOLD);
  }

  public AlphaFoldStructureParser(double contactThreshold) {
    if (contactThreshold <= 0) {
      throw new IllegalArgumentException("Contact threshold must be positive");
    }
    this.contactThreshold = contactThreshold;
  }

  public AlphaFoldStructure parse(String pdbPath) throws IOException {
    return parse(Path.of(pdbPath));
  }

  public AlphaFoldStructure parse(Path pdbPath) throws IOException {
    Objects.requireNonNull(pdbPath, "pdbPath");
    return parse(pdbPath.toFile());
  }

  public AlphaFoldStructure parse(File pdbFile) throws IOException {
    Objects.requireNonNull(pdbFile, "pdbFile");

    PDBFileReader reader = new PDBFileReader();
    FileParsingParameters params = new FileParsingParameters();
    params.setParseSecStruc(true);
    params.setLoadChemCompInfo(true);
    reader.setFileParsingParameters(params);

    Structure structure = reader.getStructure(pdbFile);

    List<Group> residueGroups = collectResidues(structure);
    double[][] distanceMatrix = buildDistanceMatrix(residueGroups);
    Map<Integer, List<ResidueContact>> contactMap = buildContactMap(residueGroups, distanceMatrix);
    Map<Integer, ResidueStructuralInfo> residueInfo = buildResidueInfo(residueGroups, contactMap);

    return new AlphaFoldStructure(residueInfo, distanceMatrix, contactThreshold);
  }

  private List<Group> collectResidues(Structure structure) {
    List<Group> residues = new ArrayList<>();
    for (Chain chain : structure.getChains()) {
      for (Group group : chain.getAtomGroups()) {
        if (group.getType().equals(GroupType.AMINOACID)) {
          residues.add(group);
        }
      }
    }
    return residues;
  }

  private double[][] buildDistanceMatrix(List<Group> residues) {
    int size = residues.size();
    double[][] distances = new double[size][size];

    for (int i = 0; i < size; i++) {
      Atom atomI = getReferenceAtom(residues.get(i));
      for (int j = i + 1; j < size; j++) {
        Atom atomJ = getReferenceAtom(residues.get(j));
        double distance = Double.NaN;
        if (atomI != null && atomJ != null) {
          try {
            distance = Calc.getDistance(atomI, atomJ);
          } catch (StructureException e) {
            distance = Double.NaN;
          }
        }
        distances[i][j] = distance;
        distances[j][i] = distance;
      }
    }
    return distances;
  }

  private Map<Integer, List<ResidueContact>> buildContactMap(
      List<Group> residues, double[][] distanceMatrix) {
    Map<Integer, List<ResidueContact>> contacts = new HashMap<>();

    for (int i = 0; i < residues.size(); i++) {
      int indexI = residues.get(i).getResidueNumber().getSeqNum();
      for (int j = i + 1; j < residues.size(); j++) {
        int indexJ = residues.get(j).getResidueNumber().getSeqNum();
        double distance = distanceMatrix[i][j];
        if (!Double.isNaN(distance) && distance <= contactThreshold) {
          contacts.computeIfAbsent(indexI, ignored -> new ArrayList<>())
              .add(new ResidueContact(indexJ, distance));
          contacts.computeIfAbsent(indexJ, ignored -> new ArrayList<>())
              .add(new ResidueContact(indexI, distance));
        }
      }
    }

    return contacts;
  }

  private Map<Integer, ResidueStructuralInfo> buildResidueInfo(
      List<Group> residues, Map<Integer, List<ResidueContact>> contactMap) {
    Map<Integer, ResidueStructuralInfo> info = new HashMap<>();

    for (Group group : residues) {
      int residueIndex = group.getResidueNumber().getSeqNum();
      String residueName = group.getPDBName();
      double confidence = extractConfidence(group);
      SecondaryStructureType secondaryStructure = extractSecondaryStructure(group);
      List<ResidueContact> contacts = contactMap.getOrDefault(residueIndex, List.of());

      info.put(
          residueIndex,
          new ResidueStructuralInfo(
              residueIndex,
              residueName,
              secondaryStructure,
              confidence,
              contacts));
    }

    return info;
  }

  private SecondaryStructureType extractSecondaryStructure(Group group) {
    Object property = group.getProperty(Group.SEC_STRUC);
    if (property instanceof SecStrucInfo secStrucInfo) {
      SecStrucType type = secStrucInfo.getType();
      if (type != null) {
        return SecondaryStructureType.fromCode(type.getTypeChar());
      }
    }
    return SecondaryStructureType.COIL;
  }

  private double extractConfidence(Group group) {
    Atom atom = getReferenceAtom(group);
    return atom != null ? atom.getTempFactor() : Double.NaN;
  }

  private Atom getReferenceAtom(Group group) {
    if (group.hasAtom("CA")) {
      return group.getAtom("CA");
    }
    if (!group.getAtoms().isEmpty()) {
      return group.getAtoms().get(0);
    }
    return null;
  }
}
