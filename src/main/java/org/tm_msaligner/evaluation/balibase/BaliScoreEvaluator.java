package org.tm_msaligner.evaluation.balibase;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Java reimplementation of the BAliBASE bali_score SP/TC evaluation logic.
 *
 * <p>The evaluator maps every non-gap residue in the test alignment back to its reference
 * alignment column and then computes:</p>
 * <ul>
 *   <li>SP score: fraction of reference residue pairs recovered in the test alignment.</li>
 *   <li>TC score: fraction of evaluated test columns that exactly recover a reference column.</li>
 * </ul>
 *
 * <p>Reference XML files use the BAliBASE {@code coreblock} column-score when present;
 * otherwise columns with at least 20% gaps are excluded, matching the original C program.</p>
 */
public final class BaliScoreEvaluator {
    private static final double DEFAULT_GAP_CUTOFF_FRACTION = 0.20;

    public BaliScoreResult evaluate(Path referenceAlignment, Path testAlignment) throws IOException {
        Alignment reference = AlignmentParsers.readReference(referenceAlignment);
        Alignment test = AlignmentParsers.readTest(testAlignment);
        return evaluate(reference, test);
    }

    public BaliScoreResult evaluate(Alignment reference, Alignment test) {
        validateSequenceSets(reference, test);

        int[] testToReference = buildTestToReferenceOrder(reference, test);
        int[] eligibleResiduesByReferenceColumn = eligibleReferenceColumns(reference);
        int maxSpPairs = calculateMaxSpPairs(eligibleResiduesByReferenceColumn);
        if (maxSpPairs <= 0) {
            throw new IllegalArgumentException("Reference alignment has no evaluable residue pairs");
        }

        Map<String, int[]> referenceResidueCodes = buildReferenceResidueCodes(reference, eligibleResiduesByReferenceColumn);
        int[][] seqCode = mapTestResiduesToReferenceColumns(reference, test, testToReference, referenceResidueCodes,
                eligibleResiduesByReferenceColumn);

        ScoreAccumulator accumulator = scoreColumns(seqCode, eligibleResiduesByReferenceColumn, test.maxLength(), reference.sequenceCount());
        double spScore = accumulator.rawSpPairs / (double) maxSpPairs;
        double tcScore = accumulator.evaluatedColumns == 0 ? 0.0 : accumulator.correctColumns / (double) accumulator.evaluatedColumns;
        return new BaliScoreResult(spScore, tcScore, accumulator.rawSpPairs, maxSpPairs,
                accumulator.correctColumns, accumulator.evaluatedColumns);
    }

    private void validateSequenceSets(Alignment reference, Alignment test) {
        if (reference.sequenceCount() != test.sequenceCount()) {
            throw new IllegalArgumentException("Different number of sequences: reference=" + reference.sequenceCount()
                    + ", test=" + test.sequenceCount());
        }
        Map<String, AlignmentSequence> testByName = new HashMap<>();
        for (AlignmentSequence sequence : test.sequences()) {
            testByName.put(sequence.comparableName(), sequence);
        }
        for (AlignmentSequence sequence : reference.sequences()) {
            AlignmentSequence testSequence = testByName.get(sequence.comparableName());
            if (testSequence == null) {
                throw new IllegalArgumentException("Sequence " + sequence.name() + " was not found in test alignment");
            }
            if (sequence.residueCount() != testSequence.residueCount()) {
                throw new IllegalArgumentException("Sequence " + sequence.name() + " has different residue counts: reference="
                        + sequence.residueCount() + ", test=" + testSequence.residueCount());
            }
        }
    }

    private int[] buildTestToReferenceOrder(Alignment reference, Alignment test) {
        Map<String, Integer> referenceIndexByName = new HashMap<>();
        for (int i = 0; i < reference.sequenceCount(); i++) {
            referenceIndexByName.put(reference.sequence(i).comparableName(), i);
        }
        int[] mapping = new int[test.sequenceCount()];
        for (int i = 0; i < test.sequenceCount(); i++) {
            Integer referenceIndex = referenceIndexByName.get(test.sequence(i).comparableName());
            if (referenceIndex == null) {
                throw new IllegalArgumentException("Sequence " + test.sequence(i).name() + " was not found in reference alignment");
            }
            mapping[i] = referenceIndex;
        }
        return mapping;
    }

    private int[] eligibleReferenceColumns(Alignment reference) {
        int[] eligible = new int[reference.maxLength()];
        int[] coreBlockMask = reference.coreBlockMask();
        if (coreBlockMask != null && coreBlockMask.length > 0) {
            int limit = Math.min(reference.maxLength(), coreBlockMask.length);
            for (int column = 0; column < limit; column++) {
                eligible[column] = coreBlockMask[column] == 1 ? reference.sequenceCount() : 0;
            }
            return eligible;
        }

        int cutoff = (int) (reference.sequenceCount() * DEFAULT_GAP_CUTOFF_FRACTION);
        if (cutoff < 1) {
            cutoff = 1;
        }
        for (int column = 0; column < reference.maxLength(); column++) {
            int gaps = 0;
            for (AlignmentSequence sequence : reference.sequences()) {
                if (sequence.isGapAt(column)) {
                    gaps++;
                }
            }
            eligible[column] = gaps >= cutoff ? 0 : reference.sequenceCount() - gaps;
        }
        return eligible;
    }

    private int calculateMaxSpPairs(int[] eligibleResiduesByReferenceColumn) {
        int max = 0;
        for (int residues : eligibleResiduesByReferenceColumn) {
            if (residues > 1) {
                max += residues * (residues - 1) / 2;
            }
        }
        return max;
    }

    private Map<String, int[]> buildReferenceResidueCodes(Alignment reference, int[] eligibleResiduesByReferenceColumn) {
        Map<String, int[]> codesByName = new HashMap<>();
        for (AlignmentSequence sequence : reference.sequences()) {
            int residueIndex = 0;
            int[] codes = new int[sequence.residueCount()];
            for (int column = 0; column < reference.maxLength(); column++) {
                if (!sequence.isGapAt(column)) {
                    codes[residueIndex++] = eligibleResiduesByReferenceColumn[column] > 0 ? column + 1 : 0;
                }
            }
            codesByName.put(sequence.comparableName(), codes);
        }
        return codesByName;
    }

    private int[][] mapTestResiduesToReferenceColumns(Alignment reference, Alignment test, int[] testToReference,
                                                       Map<String, int[]> referenceResidueCodes,
                                                       int[] eligibleResiduesByReferenceColumn) {
        int[][] seqCode = new int[reference.sequenceCount()][test.maxLength()];
        for (int testIndex = 0; testIndex < test.sequenceCount(); testIndex++) {
            AlignmentSequence testSequence = test.sequence(testIndex);
            int referenceIndex = testToReference[testIndex];
            int[] residueCodes = referenceResidueCodes.get(testSequence.comparableName());
            int residueIndex = 0;
            for (int column = 0; column < test.maxLength(); column++) {
                if (testSequence.isGapAt(column)) {
                    seqCode[referenceIndex][column] = 0;
                } else {
                    int referenceColumnCode = residueCodes[residueIndex++];
                    if (referenceColumnCode > 0 && eligibleResiduesByReferenceColumn[referenceColumnCode - 1] > 0) {
                        seqCode[referenceIndex][column] = referenceColumnCode;
                    }
                }
            }
        }
        return seqCode;
    }

    private ScoreAccumulator scoreColumns(int[][] seqCode, int[] eligibleResiduesByReferenceColumn,
                                          int testLength, int sequenceCount) {
        int rawSpPairs = 0;
        int correctColumns = 0;
        int evaluatedColumns = 0;

        for (int column = 0; column < testLength; column++) {
            Map<Integer, Integer> occurrencesByReferenceColumn = new HashMap<>();
            for (int sequence = 0; sequence < sequenceCount; sequence++) {
                int referenceColumnCode = seqCode[sequence][column];
                if (referenceColumnCode > 0) {
                    occurrencesByReferenceColumn.merge(referenceColumnCode, 1, Integer::sum);
                }
            }
            if (occurrencesByReferenceColumn.isEmpty()) {
                continue;
            }
            evaluatedColumns++;
            for (Map.Entry<Integer, Integer> entry : occurrencesByReferenceColumn.entrySet()) {
                int referenceColumnCode = entry.getKey();
                int count = Math.min(entry.getValue(), eligibleResiduesByReferenceColumn[referenceColumnCode - 1]);
                if (count > 1) {
                    rawSpPairs += count * (count - 1) / 2;
                }
            }
            if (occurrencesByReferenceColumn.size() == 1) {
                Map.Entry<Integer, Integer> only = occurrencesByReferenceColumn.entrySet().iterator().next();
                int expectedResidues = eligibleResiduesByReferenceColumn[only.getKey() - 1];
                if (only.getValue() >= expectedResidues) {
                    correctColumns++;
                }
            }
        }
        return new ScoreAccumulator(rawSpPairs, correctColumns, evaluatedColumns);
    }

    public static String formatAsBaliScoreLine(Path testAlignment, BaliScoreResult result) {
        return String.format(Locale.ROOT, "%s %.3f %.3f", testAlignment, result.spScore(), result.tcScore());
    }

    private static final class ScoreAccumulator {
        private final int rawSpPairs;
        private final int correctColumns;
        private final int evaluatedColumns;

        private ScoreAccumulator(int rawSpPairs, int correctColumns, int evaluatedColumns) {
            this.rawSpPairs = rawSpPairs;
            this.correctColumns = correctColumns;
            this.evaluatedColumns = evaluatedColumns;
        }
    }
}
