package org.tm_msaligner.score.impl;


import org.tm_msaligner.score.Score;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.util.AA;

public class LDDTStructuralScore implements Score {

  
  private final float R0;                
  private final float[] T = new float[]{0.5f, 1.0f, 2.0f, 4.0f};

  public LDDTStructuralScore(float R0) {
      this.R0 = R0;
  }

  @Override
  public <S extends TM_MSASolution> double compute(S solution, AA[][] msa) {

    int M = solution.variables().size();
    int L = solution.getAlignmentLength();

    double sumPairs = 0.0;
    int countPairs = 0;

    for (int a = 0; a < M - 1; a++) {
      float[][] DA = ((StructuralTM_MSASolution) solution).getDistanceMatrixByIndex(a);
      if (DA == null) continue;

      for (int b = a + 1; b < M; b++) {
        float[][] DB = ((StructuralTM_MSASolution) solution).getDistanceMatrixByIndex(b);
        if (DB == null) continue;

        // score del par
        double sumCols = 0.0;
        int countCols = 0;

        for (int k = 0; k < L; k++) {
          AA aaA = msa[a][k];
          AA aaB = msa[b][k];

          if (aaA.isGap() || aaB.isGap()) continue;

          int iA = aaA.getStructIndex();
          int iB = aaB.getStructIndex();
          if (iA < 0 || iB < 0) continue;

          // lDDT local para columna k (usando vecinos por columnas l)
          double colScore = lddtColumn(DA, DB, msa[a], msa[b], iA, iB, k, L);
          if (colScore >= 0.0) { // -1 si no hubo comparables
            sumCols += colScore;
            countCols++;
          }
        }

        if (countCols > 0) {
          sumPairs += (sumCols / countCols);
          countPairs++;
        }
      }
    }

    if (countPairs == 0) return 0.0;
    return sumPairs / countPairs; // en [0,1], más alto mejor
  }

  private double lddtColumn(float[][] DA, float[][] DB,
                            AA[] alnA, AA[] alnB,
                            int iA, int iB, int k, int L) {

    // counts por umbral
    int[] preserved = new int[T.length];
    int[] comparable = new int[T.length];

    for (int l = 0; l < L; l++) {
      if (l == k) continue;

      AA a2 = alnA[l];
      AA b2 = alnB[l];
      if (a2.isGap() || b2.isGap()) continue;

      int jA = a2.getStructIndex();
      int jB = b2.getStructIndex();
      if (jA < 0 || jB < 0) continue;

      float dA = DA[iA][jA];
      if (dA > R0) continue; // vecindad definida en A

      float dB = DB[iB][jB];
      float diff = Math.abs(dA - dB);

      for (int t = 0; t < T.length; t++) {
        comparable[t]++;
        if (diff <= T[t]) preserved[t]++;
      }
    }

    // si no hay comparables, no contribuye
    int validTs = 0;
    double sumFrac = 0.0;

    for (int t = 0; t < T.length; t++) {
      if (comparable[t] > 0) {
        sumFrac += (preserved[t] / (double) comparable[t]);
        validTs++;
      }
    }

    if (validTs == 0) return -1.0;
    return sumFrac / validTs;
  }

  @Override
  public boolean isAMinimizationScore() { return false; }

  @Override
  public String name() { return "lDDT_structural"; }

  @Override
  public String description() {
    return "Reference-free lDDT structural consistency over MSA" ;
  }

  public String getDescription() { return description(); }

  @Override
  public String getName() {
    return name();
  }


}
