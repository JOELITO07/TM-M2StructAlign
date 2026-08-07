package org.tm_msaligner.evaluation;

/**
 * Immutable result of comparing one multiple sequence alignment against a
 * reference alignment.
 */
public record AlignmentEvaluationResult(
    int sequenceCount,
    long referencePairCount,
    long testPairCount,
    long correctPairCount,
    double pairPrecision,
    double pairRecall,
    double pairF1,
    long referenceTmPairCount,
    long testTmPairCount,
    long correctTmPairCount,
    double tmPairPrecision,
    double tmPairRecall,
    double tmPairF1,
    long tmGapCount,
    long tmGapOpportunityCount,
    double tmGapRate) {
}
