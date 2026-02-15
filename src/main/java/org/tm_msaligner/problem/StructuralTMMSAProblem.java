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
import java.util.List;
import java.util.Map;


import org.uma.jmetal.util.errorchecking.JMetalException;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;

public class StructuralTMMSAProblem extends StandardTMMSAProblem {

  private List<float[][]> distanceMatricesByIndex;
  private List<Map<Integer, Integer>> seqPosToStructPosByIndex;

 
  public StructuralTMMSAProblem(String msaProblemFileName,
                             List<String> preComputedFiles,
                             String distanceDatasetDir) throws IOException {
    super(msaProblemFileName, preComputedFiles);

    if (distanceDatasetDir == null || distanceDatasetDir.trim().isEmpty()) {
        throw new JMetalException("Distance dataset directory is null/empty");
      }

    Path baseDir = Path.of(distanceDatasetDir);
    if (!Files.exists(baseDir) || !Files.isDirectory(baseDir)) {
        throw new JMetalException("Distance dataset directory does not exist or is not a directory: " + baseDir);
    }

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

   
}
