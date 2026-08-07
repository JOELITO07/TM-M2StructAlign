package org.tm_msaligner.evaluation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Command-line batch evaluator for execution/seed/MSA*.fasta result trees. */
public final class EvaluationMain {
  private static final String USAGE = """
      Usage:
        java -cp target/classes org.tm_msaligner.evaluation.EvaluationMain \\
          --reference <gpcrdb-reference.fasta> \\
          --topology <dataset_predicted_topologies.3line> \\
          --msa-root <execution-directory> \\
          [--output <evaluation.csv>] [--glob <MSA*.fasta>]

      The MSA root is searched recursively, so layouts such as
      execution/seed/MSA1.fasta are supported directly.
      """;

  private EvaluationMain() {
  }

  public static void main(String[] args) throws IOException {
    Configuration configuration = parseArguments(args);
    List<Path> msaFiles = discoverAlignments(configuration.msaRoot(), configuration.glob());
    if (msaFiles.isEmpty()) {
      throw new IOException(
          "No files matching '" + configuration.glob() + "' were found under "
              + configuration.msaRoot());
    }
    validateOutputTarget(configuration, msaFiles);

    AlignmentEvaluator evaluator = new AlignmentEvaluator(
        configuration.reference(), configuration.topology());
    Path outputParent = configuration.output().toAbsolutePath().getParent();
    if (outputParent != null) {
      Files.createDirectories(outputParent);
    }

    try (BufferedWriter writer = Files.newBufferedWriter(
        configuration.output(), StandardCharsets.UTF_8)) {
      writer.write(csvHeader());
      writer.newLine();

      for (Path msaFile : msaFiles) {
        AlignmentEvaluationResult result = evaluator.evaluate(msaFile);
        writer.write(toCsvRow(configuration.msaRoot(), msaFile, result));
        writer.newLine();
        System.out.printf(
            Locale.ROOT,
            "%s  Pair-F1=%.6f  TM-Pair-F1=%.6f  TM-gap-rate=%.6f%n",
            configuration.msaRoot().relativize(msaFile),
            result.pairF1(),
            result.tmPairF1(),
            result.tmGapRate());
      }
    }

    System.out.println("Evaluated alignments: " + msaFiles.size());
    System.out.println("CSV report: " + configuration.output().toAbsolutePath());
  }

  private static void validateOutputTarget(
      Configuration configuration,
      List<Path> msaFiles) throws IOException {

    Path output = configuration.output().toAbsolutePath().normalize();
    if (output.equals(configuration.reference().toAbsolutePath().normalize())
        || output.equals(configuration.topology().toAbsolutePath().normalize())) {
      throw new IOException("The CSV output must not overwrite an input file: " + output);
    }
    for (Path msaFile : msaFiles) {
      if (output.equals(msaFile.toAbsolutePath().normalize())) {
        throw new IOException("The CSV output must not overwrite an MSA file: " + output);
      }
    }
  }

  static Configuration parseArguments(String[] args) {
    if (args == null) {
      throw new IllegalArgumentException(USAGE);
    }

    Map<String, String> values = new HashMap<>();
    for (int index = 0; index < args.length; index += 2) {
      if (index + 1 >= args.length || !args[index].startsWith("--")) {
        throw new IllegalArgumentException("Invalid arguments.\n" + USAGE);
      }
      if (values.put(args[index], args[index + 1]) != null) {
        throw new IllegalArgumentException("Duplicated option " + args[index] + ".\n" + USAGE);
      }
    }

    Path reference = requiredPath(values, "--reference");
    Path topology = requiredPath(values, "--topology");
    Path msaRoot = requiredPath(values, "--msa-root");
    Path output = Path.of(values.getOrDefault(
        "--output", msaRoot.resolve("evaluation_metrics.csv").toString()));
    String glob = values.getOrDefault("--glob", "MSA*.fasta");

    for (String option : values.keySet()) {
      if (!List.of("--reference", "--topology", "--msa-root", "--output", "--glob")
          .contains(option)) {
        throw new IllegalArgumentException("Unknown option " + option + ".\n" + USAGE);
      }
    }

    if (!Files.isRegularFile(reference)) {
      throw new IllegalArgumentException("Reference file does not exist: " + reference);
    }
    if (!Files.isRegularFile(topology)) {
      throw new IllegalArgumentException("Topology file does not exist: " + topology);
    }
    if (!Files.isDirectory(msaRoot)) {
      throw new IllegalArgumentException("MSA root is not a directory: " + msaRoot);
    }
    if (glob.isBlank()) {
      throw new IllegalArgumentException("--glob must not be blank");
    }

    return new Configuration(reference, topology, msaRoot, output, glob);
  }

  private static Path requiredPath(Map<String, String> values, String option) {
    String value = values.get(option);
    if (value == null || value.isBlank()) {
      throw new IllegalArgumentException("Missing required option " + option + ".\n" + USAGE);
    }
    return Path.of(value);
  }

  static List<Path> discoverAlignments(Path msaRoot, String glob) throws IOException {
    PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:" + glob);
    List<Path> matches = new ArrayList<>();
    try (var paths = Files.walk(msaRoot)) {
      paths.filter(Files::isRegularFile)
          .filter(path -> matcher.matches(path.getFileName()))
          .forEach(matches::add);
    }
    matches.sort(Comparator.comparing(path -> msaRoot.relativize(path).toString()));
    return matches;
  }

  private static String csvHeader() {
    return String.join(",",
        "relative_path",
        "run_id",
        "msa_file",
        "sequence_count",
        "pair_precision",
        "pair_recall",
        "pair_f1",
        "tm_pair_precision",
        "tm_pair_recall",
        "tm_pair_f1",
        "tm_gap_rate",
        "reference_pairs",
        "test_pairs",
        "correct_pairs",
        "reference_tm_pairs",
        "test_tm_pairs",
        "correct_tm_pairs",
        "tm_gap_count",
        "tm_gap_opportunities");
  }

  private static String toCsvRow(
      Path msaRoot,
      Path msaFile,
      AlignmentEvaluationResult result) {

    Path relative = msaRoot.relativize(msaFile);
    String runId = relative.getNameCount() > 1
        ? relative.getParent().getFileName().toString()
        : "";

    return String.join(",",
        csv(relative.toString()),
        csv(runId),
        csv(msaFile.getFileName().toString()),
        Integer.toString(result.sequenceCount()),
        decimal(result.pairPrecision()),
        decimal(result.pairRecall()),
        decimal(result.pairF1()),
        decimal(result.tmPairPrecision()),
        decimal(result.tmPairRecall()),
        decimal(result.tmPairF1()),
        decimal(result.tmGapRate()),
        Long.toString(result.referencePairCount()),
        Long.toString(result.testPairCount()),
        Long.toString(result.correctPairCount()),
        Long.toString(result.referenceTmPairCount()),
        Long.toString(result.testTmPairCount()),
        Long.toString(result.correctTmPairCount()),
        Long.toString(result.tmGapCount()),
        Long.toString(result.tmGapOpportunityCount()));
  }

  private static String decimal(double value) {
    return Double.isNaN(value) ? "" : String.format(Locale.ROOT, "%.8f", value);
  }

  private static String csv(String value) {
    return '"' + value.replace("\"", "\"\"") + '"';
  }

  record Configuration(
      Path reference,
      Path topology,
      Path msaRoot,
      Path output,
      String glob) {
  }
}
