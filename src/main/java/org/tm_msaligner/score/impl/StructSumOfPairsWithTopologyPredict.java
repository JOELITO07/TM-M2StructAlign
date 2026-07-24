package org.tm_msaligner.score.impl;

import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.substitutionmatrix.impl.Blosum62;
import org.tm_msaligner.util.substitutionmatrix.impl.Phat;

/**
 * Topology-aware sum-of-pairs score for structural MSA solutions.
 *
 * <p>Gap scoring is symmetric with respect to the sequence order. Gap-gap
 * columns are ignored in the pairwise projection, and independent affine-gap
 * states are maintained for gaps in sequence A and sequence B.</p>
 */
public class StructSumOfPairsWithTopologyPredict implements StructuralScore {

  private enum GapState {
    NONE,
    GAP_IN_A,
    GAP_IN_B
  }

  public Phat phatMatrix;
  public Blosum62 blosum62Matrix;
  public double weightGapOpenTM;
  public double weightGapOpenNonTM;
  public double weightGapExtendTM;
  public double weightGapExtendNonTM;

  public StructSumOfPairsWithTopologyPredict(
      Phat phatMatrix,
      Blosum62 blosum62Matrix,
      double weightGapOpenTM,
      double weightGapExtendTM,
      double weightGapOpenNonTM,
      double weightGapExtendNonTM) {
    this.phatMatrix = phatMatrix;
    this.blosum62Matrix = blosum62Matrix;
    this.weightGapOpenTM = weightGapOpenTM;
    this.weightGapExtendTM = weightGapExtendTM;
    this.weightGapOpenNonTM = weightGapOpenNonTM;
    this.weightGapExtendNonTM = weightGapExtendNonTM;
  }

  @Override
  public <S extends StructuralTM_MSASolution> double compute(
      S solution, AA[][] decodedSequences) {

    int alignmentLength = solution.getAlignmentLength();
    int numberOfSequences = solution.variables().size();
    double score = 0.0;

    for (int i = 0; i < numberOfSequences - 1; i++) {
      for (int j = i + 1; j < numberOfSequences; j++) {
        GapState previousGapState = GapState.NONE;

        for (int column = 0; column < alignmentLength; column++) {
          AA aaA = decodedSequences[i][column];
          AA aaB = decodedSequences[j][column];
          boolean gapA = aaA.isGap();
          boolean gapB = aaB.isGap();

          // Gap-gap columns do not exist in the projected pairwise alignment.
          // They neither receive a penalty nor close an already open gap.
          if (gapA && gapB) {
            continue;
          }

          if (!gapA && !gapB) {
            if (isTM(aaA) && isTM(aaB)) {
              score += phatMatrix.getDistance(aaA.getLetter(), aaB.getLetter());
            } else {
              score += blosum62Matrix.getDistance(aaA.getLetter(), aaB.getLetter());
            }
            previousGapState = GapState.NONE;
            continue;
          }

          GapState currentGapState = gapA ? GapState.GAP_IN_A : GapState.GAP_IN_B;
          boolean involvesTM = isTM(aaA) || isTM(aaB);
          boolean extension = previousGapState == currentGapState;

          if (extension) {
            score -= involvesTM ? weightGapExtendTM : weightGapExtendNonTM;
          } else {
            score -= involvesTM ? weightGapOpenTM : weightGapOpenNonTM;
          }

          previousGapState = currentGapState;
        }
      }
    }

    return score;
  }

  private boolean isTM(AA aa) {
    return aa != null && aa.getType() != null && aa.getType().isTMRegion();
  }

  @Override
  public boolean isAMinimizationScore() {
    return false;
  }

  @Override
  public String name() {
    return "Sum of pairs with topology predict";
  }

  @Override
  public String description() {
    return "Symmetric topology-aware sum of pairs with affine gap penalties";
  }

  @Override
  public String getName() {
    return "SOPwTP";
  }

  public String getDescription() {
    return description();
  }
}
