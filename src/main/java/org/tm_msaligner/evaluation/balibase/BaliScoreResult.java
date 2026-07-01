package org.tm_msaligner.evaluation.balibase;

import java.util.Locale;

/** Result returned by the Java BAliBASE-compatible SP/TC evaluator. */
public final class BaliScoreResult {
    private final double spScore;
    private final double tcScore;
    private final int rawSpPairs;
    private final int maxSpPairs;
    private final int correctColumns;
    private final int evaluatedColumns;

    public BaliScoreResult(double spScore, double tcScore, int rawSpPairs, int maxSpPairs,
                           int correctColumns, int evaluatedColumns) {
        this.spScore = spScore;
        this.tcScore = tcScore;
        this.rawSpPairs = rawSpPairs;
        this.maxSpPairs = maxSpPairs;
        this.correctColumns = correctColumns;
        this.evaluatedColumns = evaluatedColumns;
    }

    public double spScore() {
        return spScore;
    }

    public double tcScore() {
        return tcScore;
    }

    public int rawSpPairs() {
        return rawSpPairs;
    }

    public int maxSpPairs() {
        return maxSpPairs;
    }

    public int correctColumns() {
        return correctColumns;
    }

    public int evaluatedColumns() {
        return evaluatedColumns;
    }

    @Override
    public String toString() {
        return String.format(Locale.ROOT,
                "SP=%.3f TC=%.3f rawSP=%d maxSP=%d correctColumns=%d evaluatedColumns=%d",
                spScore, tcScore, rawSpPairs, maxSpPairs, correctColumns, evaluatedColumns);
    }
}
