package org.tm_msaligner.mutation;

import java.util.List;

import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.util.errorchecking.Check;
import org.uma.jmetal.util.errorchecking.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/**
 * Structure-Aware Mutation: Bio-ShiftClosedGaps (rejection sampling)
 *
 * Moves an existing closed-gap block to a new position k proposed uniformly,
 * but ACCEPTS using rejection sampling:
 *   accept if u < Pacc(rk)
 * where rk is the physical residue position (1-based) mapped from aligned column k
 * for the selected sequence. If k falls on a gap for that sequence, rk = -1.
 *
 * Uses only:
 *  - gapsGroups (variables())
 *  - originalSequences (no gaps)
 *
 * Requires in TM_MSASolution:
 *  - boolean isResidueInTM(int seqIdx, int seqPos1Based)  (based on BaseType.isTMRegion())
 */
public class BioShiftClosedGapsMSAMutation implements MutationOperator<StructuralTM_MSASolution> {

  private final double mutationProbability;
  private final JMetalRandom randomGenerator;
  private final int maxTries;
  private final boolean fallbackToUniform;

  public BioShiftClosedGapsMSAMutation(double mutationProbability) {
    this(mutationProbability, 60, true);
  }

  public BioShiftClosedGapsMSAMutation(double mutationProbability, int maxTries, boolean fallbackToUniform) {
    Check.probabilityIsValid(mutationProbability);
    this.mutationProbability = mutationProbability;
    this.maxTries = Math.max(1, maxTries);
    this.fallbackToUniform = fallbackToUniform;
    this.randomGenerator = JMetalRandom.getInstance();
  }

  @Override
  public StructuralTM_MSASolution execute(StructuralTM_MSASolution solution) {
    if (solution == null) throw new JMetalException("Null parameter");
    doMutation(solution);
    return solution;
  }

  public double mutationProbability() {
    return mutationProbability;
  }

  public void doMutation(StructuralTM_MSASolution solution) {
    if (randomGenerator.nextDouble() >= mutationProbability) return;

    int selectedSequence = getSelectedSequenceWithClosedGaps(solution);

    int originalSizeAlignment = solution.getAlignmentLength(); // global alignment length
    List<Integer> gapsGroup = solution.variables().get(selectedSequence);

    int posClosedGaps = getClosedGaps(gapsGroup);

    Integer start = gapsGroup.get(posClosedGaps);
    Integer numberOfGaps = gapsGroup.get(posClosedGaps + 1) - start + 1;

    // Remove selected gap block
    gapsGroup.remove(posClosedGaps);
    gapsGroup.remove(posClosedGaps);

    // Shift remaining groups left
    for (int i = posClosedGaps; i < gapsGroup.size() - 1; i += 2) {
      gapsGroup.set(i, gapsGroup.get(i) - numberOfGaps);
      gapsGroup.set(i + 1, gapsGroup.get(i + 1) - numberOfGaps);
    }

    int maxNewPos = originalSizeAlignment - numberOfGaps;
    if (maxNewPos < 0) {
      // abort safely
      solution.variables().set(selectedSequence, gapsGroup);
      return;
    }

    Integer newpos = proposeNewPosStructureAware(solution, selectedSequence, gapsGroup, start, maxNewPos);

    // Insert the gap block at newpos (same logic as original)
    boolean added = false;
    for (int i = 0; i < gapsGroup.size() - 1; i += 2) {

      if (newpos >= gapsGroup.get(i) && newpos <= gapsGroup.get(i + 1)) {
        // extend existing gap group
        gapsGroup.set(i + 1, gapsGroup.get(i + 1) + numberOfGaps);
        added = true;

      } else if (gapsGroup.get(i) > newpos) {
        // shift subsequent groups right
        gapsGroup.set(i, gapsGroup.get(i) + numberOfGaps);
        gapsGroup.set(i + 1, gapsGroup.get(i + 1) + numberOfGaps);

        if (!added) {
          gapsGroup.add(i, newpos);
          gapsGroup.add(i + 1, newpos + numberOfGaps - 1);
          added = true;
          i += 2;
        }
      }
    }

    if (!added) {
      gapsGroup.add(newpos);
      gapsGroup.add(newpos + numberOfGaps - 1);
    }

    solution.variables().set(selectedSequence, gapsGroup);
    solution.mergeGapsGroups();
  }

  private Integer proposeNewPosStructureAware(
      StructuralTM_MSASolution solution,
      int seqIdx,
      List<Integer> gapsGroup,
      int forbiddenStart,
      int maxNewPos) {

    AAArray original = solution.getOriginalSequences().get(seqIdx);
    int originalLen = original.getSize();

    Integer accepted = null;

    for (int t = 0; t < maxTries; t++) {
      int kAligned = randomGenerator.nextInt(0, maxNewPos);
      if (kAligned == forbiddenStart) continue;

      // map aligned column -> seqPos (1-based) in UnAlignedSeq or OriginalSeq OR -1 if gap at that column
      int rk1Based = alignedPosToSeqPos1Based(kAligned, gapsGroup, originalLen);

      double pAcc = acceptanceProbability(solution, seqIdx, rk1Based, original);
      if (randomGenerator.nextDouble() < pAcc) {
        accepted = kAligned;
        break;
      }
    }

    if (accepted != null) return accepted;

    if (!fallbackToUniform) return forbiddenStart;

    Integer newpos;
    do {
      newpos = randomGenerator.nextInt(0, maxNewPos);
    } while (newpos == forbiddenStart);

    return newpos;
  }

  /**
   * Map aligned column (0-based) -> seqPos (1-based) using only gapsGroup.
   * Returns -1 if alignedPos is inside a gap block.
   */
  private int alignedPosToSeqPos1Based(int alignedPos0, List<Integer> gapsGroup, int originalLen) {
    int gapsBefore = 0;

    for (int g = 0; g < gapsGroup.size() - 1; g += 2) {
      int a = gapsGroup.get(g);
      int b = gapsGroup.get(g + 1);

      if (alignedPos0 < a) break;
      if (alignedPos0 <= b) return -1;

      gapsBefore += (b - a + 1);
    }

    int seqIndex0 = alignedPos0 - gapsBefore;
    if (seqIndex0 < 0 || seqIndex0 >= originalLen) return -1; // robust safety
    return seqIndex0 + 1; // 1-based
  }

  /**
   * Pacc(rk): prefer loops / missing-structure residues, avoid TM.
   */
  private double acceptanceProbability(StructuralTM_MSASolution solution, int seqIdx, int rk1Based, AAArray original) {

    // If the target column is already a GAP for this sequence, allow moderate acceptance
    // (placing a gap block into a column that is already a gap is usually less disruptive).
    if (rk1Based < 0) return 0.60;

    AA aa = original.AAAt(rk1Based - 1);

    // No structural representation => flexible/disordered => high acceptance
    if (aa.getStructIndex() < 0) return 0.95;

    // TM residue => strongly avoid
    if (solution.isResidueInTM(seqIdx, rk1Based)) return 0.05;

    // Non-TM residue => prefer
    return 0.90;
  }

  // ---------------- Original helpers (unchanged) ----------------

  public Integer getClosedGaps(List<Integer> gapsGroup) {
    Integer posClosedGaps = randomGenerator.nextInt(0, (gapsGroup.size() / 2) - 1);
    posClosedGaps *= 2;
    return posClosedGaps;
  }

  public int getSelectedSequenceWithClosedGaps(StructuralTM_MSASolution solution) {
    int selectedSequence;
    do {
      selectedSequence = randomGenerator.nextInt(0, solution.variables().size() - 1);
    } while (solution.variables().get(selectedSequence).isEmpty()
        || !hasSequenceClosedGaps(solution.variables().get(selectedSequence)));
    return selectedSequence;
  }

  public boolean hasSequenceClosedGaps(List<Integer> gapsGroup) {
    for (int i = 0; i < gapsGroup.size() - 1; i += 2) {
      if ((gapsGroup.get(i + 1) - gapsGroup.get(i)) > 0) return true;
    }
    return false;
  }
}
