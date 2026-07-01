package org.tm_msaligner.evaluation.balibase;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/** Immutable multiple alignment container used by the BAliBASE evaluator. */
public final class Alignment {
    private final List<AlignmentSequence> sequences;
    private final int maxLength;
    private final int[] coreBlockMask;

    public Alignment(List<AlignmentSequence> sequences, int[] coreBlockMask) {
        if (sequences == null || sequences.isEmpty()) {
            throw new IllegalArgumentException("An alignment must contain at least one sequence");
        }
        List<AlignmentSequence> copy = new ArrayList<>(sequences.size());
        int max = 0;
        for (AlignmentSequence sequence : sequences) {
            AlignmentSequence checked = Objects.requireNonNull(sequence, "sequence");
            copy.add(checked);
            max = Math.max(max, checked.alignedSequence().length());
        }
        this.sequences = Collections.unmodifiableList(copy);
        this.maxLength = max;
        this.coreBlockMask = coreBlockMask == null ? null : coreBlockMask.clone();
    }

    public List<AlignmentSequence> sequences() {
        return sequences;
    }

    public int sequenceCount() {
        return sequences.size();
    }

    public int maxLength() {
        return maxLength;
    }

    public int[] coreBlockMask() {
        return coreBlockMask == null ? null : coreBlockMask.clone();
    }

    public AlignmentSequence sequence(int index) {
        return sequences.get(index);
    }
}
