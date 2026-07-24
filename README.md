# Multiobjective Optimization for Structural-Aware Sequence Alignment of Transmembrane Proteins

_An optimization framework integrating AlphaFold2-guided structural constraints._

## 📌 Overview

This project implements a transmembrane protein sequence alignment framework combining:

- **Sequence information** (similarities, clustering).
- **Structural data** derived from AlphaFold2 predictions and topology comparisons.
- **Multiobjective optimization** to balance accuracy, stability, and complexity.

The major innovation is the **Structural Aware Alignment** module, which incorporates structural data into alignments and evaluates solutions with metrics such as TM-score, RMSD, or predicted topologies.

## 🔧 Key Features

1. **Alignment strategies**
   - Reference alignments (`.msf`, `.fasta`) and support for external tools (MAFFT, Kalign, ClustalW, etc.).
   - Custom `TM-M2StructAlign` module with mutation and crossover operators tailored to transmembrane sequences.

2. **Structural Aware Alignment (SAA)**
   - Loading predicted topology files (`*_predicted_topologies.3line`).
   - Penalties/bonuses based on agreement between topology and generated structure.
   - Utilization of AlphaFold2 data: model superposition and TM domain extraction.
   - Specialized objectives: topology concordance, helix overlap, Cα distance.

3. **Multiobjective optimization**
   - Evolutionary algorithms run over populations of MSAs.
   - Pareto front evaluation and non-dominated individuals.
   - Metrics included: Baliscore, gap weight, structural similarity.

4. **Benchmarking and testing**
   - Datasets in `resources/benchmarks/ref7` and `custom_tests`.
   - Precomputed results in `precomputed_solutions`.
   - Test scripts under `tests/` for integrity and method comparison.

5. **Lightweight web interface**
   - MSA browser with `libs/msabrowser.js` and `style.css`.
   - Interactive visualization and exploration of results.

6. **Result generation and analysis**
   - Outputs in TSV, FASTA, and MSF formats.
   - Computation of reference fronts (`referenceFronts/*.csv`).

## 🛠 Installation

```bash
# Requires Java 17+ and Maven
git clone https://.../TM-M2StructAlign.git
cd TM-M2StructAlign
mvn clean package
```

The runnable JAR will be in `target/` and can be executed with:

```bash
java -jar target/tm-msaligner.jar [options]
```

## 🚀 Usage

Example executions:

```bash
# standard multiobjective alignment
java -jar tm-msaligner.jar \
  -input resources/benchmarks/ref7/7tm/7tm.tfa \
  -topology resources/benchmarks/ref7/7tm/7tm_predicted_topologies.3line \
  -mode structural-aware \
  -generations 1000

# generate TSV report
java -jar tm-msaligner.jar -report results/FUN.tsv
```

Parameters include:

- `-mode structural-aware` to enable structural constraints.
- `-topology` path to predicted topology file from AlphaFold2.
- `-obj` multiple objectives configurable (`baliscore`, `topology`, …).

Refer to the documentation in `src/main/java/org/...` for additional details.

## 📁 Repository Structure

- `src/main/java/org/…` – main source code.
- `resources/…` – example data, benchmarks, and tests.
- `precomputed_solutions/` – reference solutions for comparison.
- `tests/` – unit and integration tests.

## ✨ Structural Aware Alignment Changes

- Initial implementation of the **structural evaluation function**.
- Reading and normalization of 3-line topology format.
- Adaptation of genetic operators to consider TM information.
- Dedicated module `AlignmentStructuralEvaluator` with superposition algorithms.
- Configuration interface for selecting “hard” vs “soft” constraints.
- Added metrics: helix agreement, intra-RMSD, and TM-score.
- Performance improvements via distance caches and parallelism.

## ✅ Example Results

- `resources/benchmarks/*` contain reference alignments.
- `custom_tests/` allows reproducible experiments with custom data.
- Pareto fronts can be visualized with the web tool.

## 🧪 Adding New Data

1. Copy `*.tfa` sequences and `*_predicted_topologies.3line` into a new folder under `resources/`.
2. Add the path to the configuration file or provide it via command line.
3. Run the JAR with `structural-aware` mode.

## 📄 License & Credits

Project licensed under the [LICENSE](LICENSE).  
Developed as part of the transmembrane structural alignment initiative.

---

> 📝 **Note:** The project’s focus is on the Structural Aware Alignment module; the README details its components, activation, and recent modifications.

