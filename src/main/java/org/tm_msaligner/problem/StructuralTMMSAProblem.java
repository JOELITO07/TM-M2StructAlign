package org.tm_msaligner.problem;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.errorchecking.JMetalException;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;

public class StructuralTMMSAProblem extends AbstractStructuralTM_MSAProblem<StructuralTM_MSASolution> {

  public List<AAArray> originalSequences;
  public List<List<AAArray>> listOfPrecomputedStringAlignments;
  public List<StringBuilder> listOfSequenceNames;
  public long[] MaxMinSegmentAlignScore;
  private List<float[][]> distanceMatricesByIndex;
  private List<Map<Integer, Integer>> seqPosToStructPosByIndex;

 
  public StructuralTMMSAProblem(String msaProblemFileName,
                             List<String> preComputedFiles,
                             String distanceDatasetDir) throws IOException {
    
    if (preComputedFiles.size() < 2) {
      throw new JMetalException(
          "Wrong number of Pre-computed Alignments, Minimum 2 files are required");
    }

    if (distanceDatasetDir == null || distanceDatasetDir.trim().isEmpty()) {
        throw new JMetalException("Distance dataset directory is null/empty");
      }

    Path baseDir = Path.of(distanceDatasetDir);
    if (!Files.exists(baseDir) || !Files.isDirectory(baseDir)) {
        throw new JMetalException("Distance dataset directory does not exist or is not a directory: " + baseDir);
    }

    readSequenceFromFile(msaProblemFileName);
    setNumberOfVariables(originalSequences.size());
    MaxMinSegmentAlignScore = getMaxMinScoreSegmentAlign();
    listOfPrecomputedStringAlignments = readPreComputedAlignments(preComputedFiles);
    
    loadDistanceAndIndexMatrices(distanceDatasetDir);
    annotateOriginalSequencesWithStructure();

  }




  private void annotateOriginalSequencesWithStructure() {
      if (seqPosToStructPosByIndex == null) return; 

      for (int m = 0; m < originalSequences.size(); m++) {
        AAArray seq = originalSequences.get(m);
        Map<Integer,Integer> map = seqPosToStructPosByIndex.get(m); 
        if (map == null) {
          for (int p = 0; p < seq.getSize(); p++) {
            AA aa = seq.AAAt(p);
            aa.setSeqIndex(p);        
            aa.setStructIndex(-1);
          }
          continue;
        }

        // IMPORTANT: map uses seqPos 1-based
        for (int p = 0; p < seq.getSize(); p++) {
          AA aa = seq.AAAt(p);
          aa.setSeqIndex(p);

          Integer structPos = map.get(p + 1); // <-- 1-based
          aa.setStructIndex(structPos != null ? structPos : -1); //Si no tiene posición estructural, se marca como -1
        }
      }
}


  private void loadDistanceAndIndexMatrices(String distanceDatasetDir) throws IOException {
      
      Path baseDir = Path.of(distanceDatasetDir);
      int n = originalSequences.size();
      distanceMatricesByIndex = new ArrayList<>(n);
      seqPosToStructPosByIndex = new ArrayList<>(n);

      for (int idx = 0; idx < listOfSequenceNames.size(); idx++) {
        String seqName = listOfSequenceNames.get(idx).toString().trim();

        Path csvD = baseDir.resolve(seqName + "_D.csv");
        Path csvI = baseDir.resolve(seqName + "_idx.csv");

        if (!Files.exists(csvD) || !Files.exists(csvI)) {
          System.out.println("[WARN] Missing structural files for '" + seqName + "'. storing nulls for Sequece " + seqName);

          distanceMatricesByIndex.add(null);
          seqPosToStructPosByIndex.add(null);
          continue;
        }

        float[][] D = null;
        int[] idxArr = null;
        Map<Integer, Integer> map = null;

        try {
          D = readDistanceCsv(csvD);
          idxArr = readIndexCsv(csvI);

          validateDistanceAndIdx(seqName, D, idxArr, originalSequences.get(idx));
          map = buildSeqPosToStructPos(seqName, idxArr);

          distanceMatricesByIndex.add(D);
          seqPosToStructPosByIndex.add(map);

        } catch (Exception e) {

          System.out.println("[WARN] Failed to load/validate structural data for '" + seqName
              + "'. Storing nulls. Reason: " + e.getMessage());

          distanceMatricesByIndex.add(null);
          seqPosToStructPosByIndex.add(null);
        }
      }
}


  private void validateDistanceAndIdx(String seqName, float[][] D, int[] idxArr, AAArray originalSeq) {
    if (D == null || D.length == 0 || D[0].length == 0) {
      throw new JMetalException("Empty distance matrix for '" + seqName + "'");
    }
    if (D.length != D[0].length) {
      throw new JMetalException("Distance matrix is not square for '" + seqName +
          "' (rows=" + D.length + ", cols=" + D[0].length + ")");
    }
    if (idxArr == null || idxArr.length == 0) {
      throw new JMetalException("Empty idx array for '" + seqName + "'");
    }
    if (D.length != idxArr.length) {
      throw new JMetalException("Mismatch D vs idx length for '" + seqName +
          "': D=" + D.length + " idx=" + idxArr.length);
    }

    int ungappedLen = getUngappedLength(originalSeq);

    int prev = Integer.MIN_VALUE;
    for (int i = 0; i < idxArr.length; i++) {
      int v = idxArr[i];
      if (v <= 0) {
        throw new JMetalException("Invalid idx (<=0) for '" + seqName + "' at structPos=" + i + " idx=" + v);
      }
      if (v < prev) {
        throw new JMetalException("Non-increasing idx for '" + seqName + "' at structPos=" + i +
            " prev=" + prev + " curr=" + v);
      }
      if (v > ungappedLen) {
        throw new JMetalException("idx out of range for '" + seqName + "': idx=" + v +
            " but ungappedLen=" + ungappedLen +
            ". (Tu PDB/idx no corresponde a la secuencia de entrada)");
      }
      prev = v;
    }

    assertSymmetric(D, seqName, 1e-3f);
  }

  private Map<Integer, Integer> buildSeqPosToStructPos(String seqName, int[] idxArr) {
    Map<Integer, Integer> map = new HashMap<>(idxArr.length * 2);

    for (int structPos = 0; structPos < idxArr.length; structPos++) {
      int seqPos = idxArr[structPos]; // 1-based

      if (!map.containsKey(seqPos)) {
        map.put(seqPos, structPos);
      }
    }

    if (map.isEmpty()) {
      throw new JMetalException("Empty seqPos->structPos map for '" + seqName + "'");
    }
    return map;
  }

  private float[][] readDistanceCsv(Path csvPath) throws IOException {
    List<String> lines = Files.readAllLines(csvPath, StandardCharsets.UTF_8);

    List<String> rows = new ArrayList<>();
    for (String line : lines) {
      if (line != null) {
        String t = line.trim();
        if (!t.isEmpty()) rows.add(t);
      }
    }

    if (rows.isEmpty()) {
      throw new JMetalException("Empty distance CSV: " + csvPath);
    }

    int n = -1;
    float[][] D = null;

    for (int r = 0; r < rows.size(); r++) {
      String[] parts = rows.get(r).split(",");

      if (n < 0) {
        n = parts.length;
        D = new float[n][n];
      } else if (parts.length != n) {
        throw new JMetalException("Inconsistent row length at row " + r + " in " + csvPath +
            " (expected " + n + " columns, got " + parts.length + ")");
      }

      if (r >= n) {
        throw new JMetalException("More rows than expected (not NxN) in " + csvPath +
            " (row index " + r + ", expected max " + (n - 1) + ")");
      }

      for (int c = 0; c < n; c++) {
        D[r][c] = Float.parseFloat(parts[c].trim());
      }
    }

    if (rows.size() != n) {
      throw new JMetalException("Matrix is not square in " + csvPath +
          " (rows=" + rows.size() + ", cols=" + n + ")");
    }

    return D;
  }

  private int[] readIndexCsv(Path idxPath) throws IOException {
    
    String content = Files.readString(idxPath, StandardCharsets.UTF_8).trim();
    if (content.isEmpty()) {
      throw new JMetalException("Empty idx CSV: " + idxPath);
    }

    String[] parts = content.split("[,\\s]+"); // coma o whitespace
    int[] idx = new int[parts.length];
    for (int i = 0; i < parts.length; i++) {
      idx[i] = Integer.parseInt(parts[i].trim());
    }
    return idx;
  }

  private int getUngappedLength(AAArray seq) {
    int len = 0;
    for (int i = 0; i < seq.getSize(); i++) {
      if (seq.AAAt(i).getLetter() != AA.GAP_IDENTIFIER) {
        len++;
      }
    }
    return len;
  }

  private void assertSymmetric(float[][] D, String seqName, float eps) {
    int n = D.length;
    for (int i = 0; i < n; i++) {
      if (Math.abs(D[i][i]) > eps) {
        throw new JMetalException("Non-zero diagonal in distance matrix for " + seqName + " at i=" + i);
      }
      for (int j = i + 1; j < n; j++) {
        if (Math.abs(D[i][j] - D[j][i]) > eps) {
          throw new JMetalException("Non-symmetric distance matrix for " + seqName +
              " at (" + i + "," + j + "): " + D[i][j] + " vs " + D[j][i]);
        }
      }
    }
  }

  public float[][] getDistanceMatrixByIndex(int seqIndex) {
    if (distanceMatricesByIndex == null || seqIndex < 0 || seqIndex >= distanceMatricesByIndex.size()) {
      throw new JMetalException("Distance matrices not loaded or invalid index: " + seqIndex);
    }
    return distanceMatricesByIndex.get(seqIndex);
  }


  /** seqPos (1-based, en secuencia sin gaps) -> structPos (0-based, fila/col en D) */
  public Map<Integer, Integer> getSeqPosToStructPosMapByIndex(int seqIndex) {
    if (seqPosToStructPosByIndex == null || seqIndex < 0 || seqIndex >= seqPosToStructPosByIndex.size()) {
      throw new JMetalException("Mapping not loaded or invalid index: " + seqIndex);
    }
    return seqPosToStructPosByIndex.get(seqIndex);
  }

   
  

  
   public List<List<AAArray>> readPreComputedAlignments(List<String> dataFiles) {
    List<List<AAArray>> listPreAlignments = new ArrayList<List<AAArray>>();
    for (String dataFile : dataFiles) {
      try {
        List<AAArray> seqAligned = readDataFromFastaFile(dataFile);
        StructuralTM_MSASolution sol = new StructuralTM_MSASolution(seqAligned, this);
        if (!sol.isValid()) {
          JMetalLogger.logger.warning("MSA in file " + dataFile + " is not Valid");
        } else {
          listPreAlignments.add(seqAligned);
        }
      } catch (Exception e) {
        throw new JMetalException(
            "Error reading data from fasta files " + dataFile + ". Message: " + e);
      }
    }

    if (listPreAlignments.size() < 2) {
      throw new JMetalException("More than one PreComputedAlignment is needed");
    }

    return listPreAlignments;
  }

  /**
   * Reads a precomputed alignment and restores the canonical input order by
   * matching normalized FASTA identifiers. A trailing .pdb extension is
   * removed from both input and precomputed identifiers, e.g.
   * ada2b_human.pdb -> ada2b_human.
   */
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

    /*
     * Remove an optional chain identifier.
     *
     * Examples:
     * ada2b_human.pdb:A -> ada2b_human.pdb
     * ada2b_human:A     -> ada2b_human
     */
    int chainSeparator = value.indexOf(':');
    if (chainSeparator >= 0) {
      value = value.substring(0, chainSeparator);
    }

    value = value.trim();
    if (value.toLowerCase().endsWith(".pdb")) {
      value = value.substring(0, value.length() - 4);
    }

    value = value.trim();

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

 /* public List<AAArray> readDataFromFastaFile(String dataFile)
      throws IOException {

    List<AAArray> sequenceList = new ArrayList<AAArray>();

    LinkedHashMap<String, ProteinSequence> sequences =
        FastaReaderHelper.readFastaProteinSequence(new File(dataFile));

    for (Map.Entry<String, ProteinSequence> entry : sequences.entrySet()) {
      sequenceList.add(new AAArray(entry.getValue().getSequenceAsString()));
    }

    return sequenceList;
  }*/

  public List<StringBuilder> getListOfSequenceNames() {
    return listOfSequenceNames;
  }

  public List<StringBuilder> readSeqNameFromAlignment(String dataFile)
      throws IOException, CompoundNotFoundException {

    List<StringBuilder> SeqNameList = new ArrayList<StringBuilder>();

    LinkedHashMap<String, ProteinSequence> sequences =
        FastaReaderHelper.readFastaProteinSequence(new File(dataFile));

    for (Map.Entry<String, ProteinSequence> entry : sequences.entrySet()) {
      SeqNameList.add(new StringBuilder(entry.getValue().getOriginalHeader()));
    }

    return SeqNameList;
  }

  public void printSequenceListToFasta(List<AAArray> solutionList, String fileName)
      throws Exception {
    List<ProteinSequence> proteinSequences = new ArrayList<ProteinSequence>();
    for (AAArray sequence : solutionList) {
      proteinSequences.add(new ProteinSequence(sequence.toString()));
    }

    FastaWriterHelper.writeProteinSequence(new File(fileName), proteinSequences);
  }

  public long[] getMaxMinScoreSegmentAlign() {
    long[] MaxMinScores = new long[2];
    long MaxScore = 0, MinScore = 0;
    AAArray seq;
    int numSeqs = originalSequences.size();

    for (int i = 0; i < numSeqs - 1; i++) {
      seq = originalSequences.get(i);
      for (int l = 0; l < seq.getSize(); l++) {
        MinScore += (-1) * (numSeqs - i - 1);
        if (seq.AAAt(l).getType().isTMRegion()) {
          MaxScore += 4 * (numSeqs - i - 1);
        } else {
          MaxScore += 2 * (numSeqs - i - 1);
        }
      }
    }
    MaxMinScores[0] = MaxScore;
    MaxMinScores[1] = MinScore;
    return MaxMinScores;
  }

  public StructuralTM_MSASolution evaluate(StructuralTM_MSASolution tm_msaSolution) { return null; }
  public StructuralTM_MSASolution createSolution() { return null; }

  void readSequenceFromFile(String file) {
    originalSequences = new ArrayList<AAArray>();
    listOfSequenceNames = new ArrayList<StringBuilder>();

    try {
      BufferedReader in = new BufferedReader(new FileReader(file));

      int status = 0;
      String line, regiones;
      int posSepName;

      for (line = in.readLine(); line != null; line = in.readLine()) {
        line = line.trim();
        if (status == 0) {
          if (line.length() > 0 && line.charAt(0) == '>') {
            posSepName = line.indexOf('|');
            listOfSequenceNames.add(new StringBuilder(
                posSepName > 0 ? line.substring(1, posSepName) : line.substring(1)));
          } else {
            throw new IOException("Name of Sequence must starts with '>'");
          }
          status = 1;
        } else if (status == 1) {
          regiones = in.readLine().trim();
          if (regiones == null) {
            throw new IOException("Regions of Sequence is empty");
          }
          if (regiones.length() != line.length()) {
            throw new IOException("Regions of Sequence is empty");
          }

          originalSequences.add(new AAArray(line, regiones));
          status = 0;
        }
      }

      if (originalSequences.size() != listOfSequenceNames.size()) {
        throw new IOException("Names wiht Sequences are not equals");
      }

    } catch (IOException e) {
      System.out.println("Error when reading " + file);
      e.printStackTrace();
    }
  }

  public int getSizeOfOriginalSequence(int i) {
    if (i < originalSequences.size()) {
      return originalSequences.get(i).getSize();
    } else {
      System.out.println("Error getting size of Original Sequence " + i +
          " and the Number of Sequences is " + originalSequences.size());
      return 0;
    }
  }

}
