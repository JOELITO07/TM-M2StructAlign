package org.tm_msaligner.score.impl;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.tm_msaligner.score.Score;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.util.AA;

public class StructurePreservingScore implements Score {

  private final String distanceMatrixPath;
  private final String weightMatrixPath;
  private final double contactThreshold;
  private final boolean useWeights;
  private final double paeTemperature;

  public StructurePreservingScore(String distDir, String weightDir, boolean useWeights) {
    this(distDir, weightDir, useWeights, 5.0, 8.0);
  }

  public StructurePreservingScore(
      String distDir,
      String weightDir,
      boolean useWeights,
      double paeTemperature,
      double contactThreshold) {
    if (distDir == null) {
      throw new IllegalArgumentException("Distance matrix directory must not be null");
    }
    if (useWeights && weightDir == null) {
      throw new IllegalArgumentException("Weight matrix directory must not be null when weights are used");
    }
    this.distanceMatrixPath = distDir;
    this.weightMatrixPath = weightDir;
    this.useWeights = useWeights;
    this.paeTemperature = paeTemperature;
    this.contactThreshold = contactThreshold;
  }

  @Override
  public String getName() {
    return "SPS";
  }

  @Override
  public String getDescription() {
    return "Structure-Preserving Score based on AlphaFold distance and PAE matrices";
  }

  @Override
  public <S extends TM_MSASolution> double compute(S solution, AA[][] decodedSequences) {
    int numSeqs = solution.variables().size();
    int alnLength = solution.getAlignmentLength();

    double totalPairs = 0.0;
    double conservedSum = 0.0;

    double[][][] distMatrices = new double[numSeqs][][];
    double[][][] weightMatrices = new double[numSeqs][][];

    for (int i = 0; i < numSeqs; i++) {
      String name = solution.getSequenceName(i);
      distMatrices[i] = loadMatrix(Paths.get(distanceMatrixPath, name + "_D.csv"));
      if (useWeights) {
        weightMatrices[i] = loadMatrix(Paths.get(weightMatrixPath, name + "_W.csv"));
      }
    }

    for (int p = 0; p < alnLength - 1; p++) {
      for (int q = p + 1; q < alnLength; q++) {
        double pairConservation = 0.0;
        double pairWeight = 0.0;
        int comparisons = 0;

        for (int a = 0; a < numSeqs - 1; a++) {
          for (int b = a + 1; b < numSeqs; b++) {
            int resP_A = mapResidueIndex(decodedSequences[a], p);
            int resQ_A = mapResidueIndex(decodedSequences[a], q);
            int resP_B = mapResidueIndex(decodedSequences[b], p);
            int resQ_B = mapResidueIndex(decodedSequences[b], q);

            if (resP_A < 0 || resQ_A < 0 || resP_B < 0 || resQ_B < 0) {
              continue;
            }

            if (!isValidPair(distMatrices[a], resP_A, resQ_A)
                || !isValidPair(distMatrices[b], resP_B, resQ_B)) {
              continue;
            }

            boolean contactA = distMatrices[a][resP_A][resQ_A] <= contactThreshold;
            boolean contactB = distMatrices[b][resP_B][resQ_B] <= contactThreshold;

            if (contactA == contactB) {
              pairConservation += 1.0;
            }

            double wA = 1.0;
            double wB = 1.0;
            if (useWeights) {
              if (isValidPair(weightMatrices[a], resP_A, resQ_A)) {
                wA = computeWeight(weightMatrices[a][resP_A][resQ_A]);
              }
              if (isValidPair(weightMatrices[b], resP_B, resQ_B)) {
                wB = computeWeight(weightMatrices[b][resP_B][resQ_B]);
              }
            }

            pairWeight += (wA + wB) / 2.0;
            comparisons++;
          }
        }

        if (comparisons > 0) {
          double meanConservation = pairConservation / comparisons;
          double meanWeight = pairWeight / comparisons;
          conservedSum += useWeights ? meanWeight * meanConservation : meanConservation;
          totalPairs++;
        }
      }
    }

    return totalPairs > 0 ? conservedSum / totalPairs : 0.0;
  }

  private double[][] loadMatrix(Path path) {
    Path normalized = path.toAbsolutePath().normalize();
    if (!Files.exists(normalized)) {
      throw new IllegalStateException("Matrix file not found: " + normalized);
    }

    List<double[]> rows = new ArrayList<>();
    int expectedColumns = -1;

    try (BufferedReader reader = Files.newBufferedReader(normalized)) {
      String line;
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        if (line.isEmpty()) {
          continue;
        }

        String[] tokens = line.split(",");
        if (tokens.length == 1) {
          tokens = line.split("\\s+");
        }

        if (expectedColumns == -1) {
          expectedColumns = tokens.length;
        } else if (tokens.length != expectedColumns) {
          throw new IllegalStateException(
              "Matrix " + normalized + " has inconsistent number of columns at row "
                  + (rows.size() + 1));
        }

        double[] row = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
          row[i] = Double.parseDouble(tokens[i].trim());
        }
        rows.add(row);
      }
    } catch (IOException e) {
      throw new IllegalStateException("Unable to read matrix from " + normalized, e);
    }

    double[][] matrix = new double[rows.size()][];
    for (int i = 0; i < rows.size(); i++) {
      matrix[i] = rows.get(i);
    }
    return matrix;
  }

  private boolean isValidPair(double[][] matrix, int i, int j) {
    return matrix != null
        && i >= 0
        && j >= 0
        && i < matrix.length
        && j < matrix[i].length;
  }

  private double computeWeight(double paeValue) {
    if (!useWeights) {
      return 1.0;
    }
    return Math.exp(-paeValue / paeTemperature);
  }

  private int mapResidueIndex(AA[] sequence, int alnPos) {
    if (sequence == null || alnPos < 0 || alnPos >= sequence.length) {
      return -1;
    }

    if (sequence[alnPos].isGap()) {
      return -1;
    }

    int index = 0;
    for (int i = 0; i < alnPos; i++) {
      if (!sequence[i].isGap()) {
        index++;
      }
    }
    return index;
  }

  @Override
  public boolean isAMinimizationScore() {
    return false;
  }

  @Override
  public String name() {
    return "Structure-Preserving Score";
  }

  @Override
  public String description() {
    return "Evaluates structural coherence using AlphaFold distances";
  }
}
