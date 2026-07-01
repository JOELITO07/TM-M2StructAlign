package org.tm_msaligner.evaluation.balibase;

import java.nio.file.Path;

/** Command-line entry point compatible with: bali_score_elm ref_aln test_aln [-v]. */
public final class BaliScoreCli {
    private BaliScoreCli() {
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2 || args.length > 3) {
            System.err.println("Usage: BaliScoreCli ref_aln test_aln [-v]");
            System.err.println("       ref_aln  reference alignment in BAliBASE XML, MSF/GCG, FASTA, or CLUSTAL format");
            System.err.println("       test_aln alignment to evaluate in MSF/GCG, FASTA, or CLUSTAL format");
            System.exit(1);
        }
        Path reference = Path.of(args[0]);
        Path test = Path.of(args[1]);
        BaliScoreEvaluator evaluator = new BaliScoreEvaluator();
        BaliScoreResult result = evaluator.evaluate(reference, test);
        System.out.println(BaliScoreEvaluator.formatAsBaliScoreLine(test, result));
        if (args.length == 3 && "-v".equals(args[2])) {
            System.out.println(result);
        }
    }
}
