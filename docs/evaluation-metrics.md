# MSA evaluation metrics

The package `org.tm_msaligner.evaluation` compares every generated MSA with a
GPCRdb reference alignment and reuses the existing `AAArray`, `AA`, and
`BaseType` topology model to identify transmembrane residues.

## Metrics

### Pair Precision, Recall and F1

For every MSA column, the evaluator records each aligned pair of non-gap
residues as `(protein A, original position A, protein B, original position B)`.
Original positions are ungapped, zero-based positions internally.

```
Pair precision = correct test pairs / test pairs
Pair recall    = correct test pairs / reference pairs
Pair F1        = 2 * precision * recall / (precision + recall)
```

This makes the comparison independent of sequence order, MSA length, and gap
placement.

### TM-Pair F1

TM-Pair Precision, Recall, and F1 use the same definitions but retain only
pairs for which both residues satisfy:

```java
aa.getType().isTMRegion()
```

Topology is loaded from the project's three-line
`*_predicted_topologies.3line` files.

### TM-gap rate

For every TM residue, the evaluator inspects the aligned character in every
other protein. A TM-gap event occurs when that character is a gap.

```
TM-gap rate = TM residue-to-gap events / (TM residues * (N - 1))
```

The denominator makes the result range from 0 to 1 and prevents datasets with
more proteins or more TM residues from receiving larger values merely because
of their size. Lower is better.

## Batch execution

Build with Java 17 and Maven:

```bash
mvn clean package
```

Evaluate all `MSA*.fasta` files recursively under a directory such as
`ejecuciones/semilla/MSA1.fasta`:

```bash
java -cp target/classes org.tm_msaligner.evaluation.EvaluationMain \
  --reference /data/GPCRdb/reference_alignments/classA_001.fasta \
  --topology /data/sequences/tmregions/classA_001_predicted_topologies.3line \
  --msa-root /results/ejecuciones/classA_001 \
  --output /results/evaluation_classA_001.csv \
  --glob 'MSA*.fasta'
```

The search is recursive. The output CSV contains the three requested metrics,
their Precision/Recall components, and all raw counts required to audit or
recompute each score.

Identifiers are matched case-insensitively. The reader also normalizes UniProt
`sp|ID|...`/`tr|ID|...` headers, chain suffixes, and `.pdb`, `.ent`, `.cif`, or
`.mmcif` extensions. Before scoring, it verifies that the normalized protein
sets and all ungapped sequences agree across the topology, reference, and test
files.
