# Java BAliBASE SP/TC evaluator

This module reimplements the core `bali_score_elm` behavior in Java so the project can compute BAliBASE-style SP and TC scores without compiling the original C program or installing Expat.

## Classes

- `org.tm_msaligner.evaluation.balibase.BaliScoreEvaluator`: programmatic API.
- `org.tm_msaligner.evaluation.balibase.BaliScoreCli`: command-line entry point.
- `org.tm_msaligner.evaluation.balibase.BaliScoreResult`: score object with SP, TC and raw counts.

## Supported input

Reference alignment:

- BAliBASE XML with `<seq-name>`, `<seq-data>` and optional `coreblock` column score.
- MSF/GCG.
- FASTA.
- CLUSTAL.

Test alignment:

- MSF/GCG.
- FASTA.
- CLUSTAL.
- XML is also parsed for convenience.

## Command-line usage

After building the Maven project:

```bash
mvn package
java -cp target/classes org.tm_msaligner.evaluation.balibase.BaliScoreCli reference.xml test.msf -v
```

Output follows the original BAliBASE style:

```text
test.msf 0.873 0.642
```

Where the first number is the normalized SP score and the second is the TC score.

## Java API usage

```java
BaliScoreEvaluator evaluator = new BaliScoreEvaluator();
BaliScoreResult result = evaluator.evaluate(
    Path.of("reference.xml"),
    Path.of("test.msf")
);

System.out.println(result.spScore());
System.out.println(result.tcScore());
```

## Notes

- Sequence names are matched case-insensitively, so the order in the reference and test files can differ.
- Residue counts must match for each sequence after removing gaps.
- XML references use BAliBASE `coreblock` columns when available. If no coreblock mask exists, columns with at least 20% gaps are excluded, following the original BAliBASE C logic.
