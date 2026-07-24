package org.tm_msaligner.evaluation.balibase;

import java.util.Locale;
import java.util.Objects;

/** A named aligned sequence. Gaps are normalized to '-'. */
public final class AlignmentSequence {
    private final String name;
    private final String alignedSequence;

    public AlignmentSequence(String name, String alignedSequence) {
        this.name = normalizeName(name);
        this.alignedSequence = normalizeSequence(alignedSequence);
    }

    public String name() {
        return name;
    }

    public String comparableName() {
        return name.toLowerCase(Locale.ROOT);
    }

    public String alignedSequence() {
        return alignedSequence;
    }

    public int length() {
        return alignedSequence.length();
    }

    public char charAtOrGap(int column) {
        return column < alignedSequence.length() ? alignedSequence.charAt(column) : '-';
    }

    public boolean isGapAt(int column) {
        return isGap(charAtOrGap(column));
    }

    public int residueCount() {
        int count = 0;
        for (int i = 0; i < alignedSequence.length(); i++) {
            if (!isGap(alignedSequence.charAt(i))) {
                count++;
            }
        }
        return count;
    }

    public static boolean isGap(char residue) {
        return residue == '-' || residue == '.' || residue == '~';
    }

    private static String normalizeName(String name) {
        String value = Objects.requireNonNull(name, "name").trim();
        if (value.isEmpty()) {
            throw new IllegalArgumentException("Sequence name cannot be empty");
        }
        return value;
    }

    private static String normalizeSequence(String sequence) {
        String value = Objects.requireNonNull(sequence, "alignedSequence");
        StringBuilder normalized = new StringBuilder(value.length());
        for (int i = 0; i < value.length(); i++) {
            char c = value.charAt(i);
            if (Character.isWhitespace(c) || Character.isDigit(c)) {
                continue;
            }
            if (c == '.' || c == '~') {
                normalized.append('-');
            } else if (c == '*') {
                normalized.append('X');
            } else if (Character.isLetter(c) || c == '-') {
                normalized.append(Character.toUpperCase(c));
            }
        }
        return normalized.toString();
    }
}
