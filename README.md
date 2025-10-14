# Multiobjective Optimization Applied to the Structure-Based Sequence Alignment of Transmembrane Proteins with AlphaFold2-Guided Constraints

## Structure-aware objective

The project now includes a `StructurePreservingScore` objective that biases the optimisation towards
alignments coherent with AlphaFold2 models. The score blends two complementary signals:

* **Contact preservation** – penalises the insertion of gaps between residues that are close in the
  3D structure (default threshold of 8 Å).
* **Secondary structure agreement** – rewards columns where helices or strands remain aligned across
  homologous proteins.

The score is available in `org.tm_msaligner.score.impl.StructurePreservingScore` and expects a map
between FASTA headers and parsed `AlphaFoldStructure` instances (see
`org.tm_msaligner.util.structure.AlphaFoldStructureParser`).

