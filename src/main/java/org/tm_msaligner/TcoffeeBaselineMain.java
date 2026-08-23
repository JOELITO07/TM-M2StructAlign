package org.tm_msaligner;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.score.impl.LDDTStructuralScore;
import org.tm_msaligner.score.impl.StructSumOfPairsWithTopologyPredict;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;
import org.tm_msaligner.util.AAArray;
import org.tm_msaligner.util.substitutionmatrix.impl.Blosum62;
import org.tm_msaligner.util.substitutionmatrix.impl.Phat;

/**
 * Evaluates every T-Coffee alignment in the GPCRdb precomputed directory with
 * the same two objectives used by {@link TM_MSAStructAlignerGPCRdb}.
 *
 * <p>Expected layout:</p>
 *
 * <pre>
 * GPCRdb/
 *   precomputed/classA_001/tcoffee_classA_001.fasta
 *   sequences/tmregions/classA_001_predicted_topologies.3line
 *   distances/classA_001/*_D.csv and *_idx.csv
 * </pre>
 *
 * <p>Outputs:</p>
 *
 * <ul>
 *   <li>{@code precomputed/baselines_tcoffee.csv}: consolidated report.</li>
 *   <li>{@code precomputed/<dataset>/baseline_tcoffee.csv}: input compatible
 *       with {@code pareto_analysis.py --baselines}.</li>
 * </ul>
 */
public final class TcoffeeBaselineMain {

  private static final Path DEFAULT_PRECOMPUTED_ROOT = Path.of(
      "C:\\GDrive2026\\TM-MSA\\Datasets\\GPCRdb\\precomputed");

  private static final double GAP_OPEN_TM = 8.0;
  private static final double GAP_EXTEND_TM = 3.0;
  private static final double GAP_OPEN_NON_TM = 3.0;
  private static final double GAP_EXTEND_NON_TM = 1.0;
  private static final float LDDT_RADIUS_ANGSTROM = 15.0f;

  private TcoffeeBaselineMain() {
    // Utility main class.
  }

  public static void main(String[] args) throws IOException {
    if (args.length > 3) {
      printUsage();
      System.exit(2);
    }

    Path precomputedRoot = args.length >= 1
        ? Path.of(args[0]).toAbsolutePath().normalize()
        : DEFAULT_PRECOMPUTED_ROOT.toAbsolutePath().normalize();

    Path inferredDataRoot = precomputedRoot.getParent();
    Path dataRoot = args.length >= 2
        ? Path.of(args[1]).toAbsolutePath().normalize()
        : inferredDataRoot;

    Path consolidatedOutput = args.length >= 3
        ? Path.of(args[2]).toAbsolutePath().normalize()
        : precomputedRoot.resolve("baselines_tcoffee.csv");

    validateDirectory(precomputedRoot, "precomputed root");
    validateDirectory(dataRoot, "GPCRdb data root");

    List<Path> datasetDirectories = listDatasetDirectories(precomputedRoot);
    if (datasetDirectories.isEmpty()) {
      throw new IOException("No dataset directories found in: " + precomputedRoot);
    }

    List<Result> results = new ArrayList<>();
    int successes = 0;

    System.out.println("Precomputed root: " + precomputedRoot);
    System.out.println("GPCRdb data root: " + dataRoot);
    System.out.println("Datasets found: " + datasetDirectories.size());
    System.out.println();

    for (Path datasetDirectory : datasetDirectories) {
      String dataset = datasetDirectory.getFileName().toString();
      Path expectedFasta = datasetDirectory.resolve("tcoffee_" + dataset + ".fasta");
      Path fasta = findCaseInsensitive(datasetDirectory, expectedFasta.getFileName().toString());

      if (fasta == null) {
        String message = "Missing " + expectedFasta.getFileName();
        System.out.println("[SKIP] " + dataset + ": " + message);
        results.add(Result.failure(dataset, "T-Coffee", expectedFasta, message));
        continue;
      }

      try {
        Result result = evaluateDataset(dataRoot, dataset, fasta);
        writeDatasetBaseline(datasetDirectory.resolve("baseline_tcoffee.csv"), result);
        results.add(result);
        successes++;
        System.out.printf(
            Locale.ROOT,
            "[OK]   %s: objective_1=%.10f objective_2=%.10f%n",
            dataset,
            result.objective1(),
            result.objective2());
      } catch (Exception exception) {
        String message = rootMessage(exception);
        System.err.println("[FAIL] " + dataset + ": " + message);
        results.add(Result.failure(dataset, "T-Coffee", fasta, message));
      }
    }

    writeConsolidated(consolidatedOutput, results);

    System.out.println();
    System.out.println("Successful datasets: " + successes + "/" + datasetDirectories.size());
    System.out.println("Consolidated baseline: " + consolidatedOutput);

    if (successes == 0) {
      System.exit(1);
    }
  }

  private static Result evaluateDataset(Path dataRoot, String dataset, Path fasta)
      throws IOException {

    Path topologyFile = dataRoot
        .resolve("sequences")
        .resolve("tmregions")
        .resolve(dataset + "_predicted_topologies.3line");
    Path distanceDirectory = dataRoot.resolve("distances").resolve(dataset);

    if (!Files.isRegularFile(topologyFile)) {
      throw new IOException("Missing topology file: " + topologyFile);
    }
    validateDirectory(distanceDirectory, "distance directory for " + dataset);

    List<StructuralScore> scores = createScores();

    /*
     * StructuralTMMSAProblem currently requires at least two precomputed
     * alignments for population initialization. Repeating the same path here
     * only satisfies construction; the target alignment is read once below
     * and evaluated directly, without running the evolutionary algorithm.
     */
    List<String> constructorAlignments = List.of(fasta.toString(), fasta.toString());
    MultiObjStructuralTMMSAProblem problem = new MultiObjStructuralTMMSAProblem(
        topologyFile.toString(),
        scores,
        constructorAlignments,
        distanceDirectory.toString(),
        dataset);

    List<AAArray> alignedSequences = problem.readDataFromFastaFile(fasta.toString());
    StructuralTM_MSASolution solution = new StructuralTM_MSASolution(alignedSequences, problem);
    if (!solution.isValid()) {
      throw new IllegalArgumentException("Invalid MSA: " + fasta);
    }

    // Match MultiObjStructuralTMMSAProblem.evaluate(): remove all-gap columns.
    solution.removeGapColumns();
    AA[][] decodedSequences = solution.decodeToMatrix();

    // Compute raw benefit values. These are the values restored before FUN.tsv
    // is written by TM_MSAStructAlignerGPCRdb.
    double objective1 = scores.get(0).compute(solution, decodedSequences);
    double objective2 = scores.get(1).compute(solution, decodedSequences);

    if (!Double.isFinite(objective1) || !Double.isFinite(objective2)) {
      throw new IllegalStateException(
          "Non-finite objective value: objective_1=" + objective1
              + ", objective_2=" + objective2);
    }

    return Result.success(dataset, "T-Coffee", fasta, objective1, objective2);
  }

  private static List<StructuralScore> createScores() {
    List<StructuralScore> scores = new ArrayList<>(2);
    scores.add(new StructSumOfPairsWithTopologyPredict(
        new Phat(8),
        new Blosum62(),
        GAP_OPEN_TM,
        GAP_EXTEND_TM,
        GAP_OPEN_NON_TM,
        GAP_EXTEND_NON_TM));
    scores.add(new LDDTStructuralScore(LDDT_RADIUS_ANGSTROM));
    return scores;
  }

  private static List<Path> listDatasetDirectories(Path precomputedRoot) throws IOException {
    List<Path> directories = new ArrayList<>();
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(precomputedRoot)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          directories.add(path);
        }
      }
    }
    directories.sort(Comparator.comparing(
        path -> path.getFileName().toString(),
        String.CASE_INSENSITIVE_ORDER));
    return directories;
  }

  private static Path findCaseInsensitive(Path directory, String expectedName) throws IOException {
    Path exact = directory.resolve(expectedName);
    if (Files.isRegularFile(exact)) {
      return exact;
    }
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(directory, "*.fasta")) {
      for (Path candidate : stream) {
        if (candidate.getFileName().toString().equalsIgnoreCase(expectedName)
            && Files.isRegularFile(candidate)) {
          return candidate;
        }
      }
    }
    return null;
  }

  private static void writeDatasetBaseline(Path output, Result result) throws IOException {
    String content = "method,objective_1,objective_2" + System.lineSeparator()
        + csv(result.method()) + ","
        + format(result.objective1()) + ","
        + format(result.objective2()) + System.lineSeparator();
    writeUtf8(output, content);
  }

  private static void writeConsolidated(Path output, List<Result> results) throws IOException {
    StringBuilder csv = new StringBuilder(
        "dataset,method,objective_1,objective_2,fasta,status,message")
        .append(System.lineSeparator());

    for (Result result : results) {
      csv.append(csv(result.dataset())).append(',')
          .append(csv(result.method())).append(',')
          .append(format(result.objective1())).append(',')
          .append(format(result.objective2())).append(',')
          .append(csv(result.fasta().toString())).append(',')
          .append(csv(result.status())).append(',')
          .append(csv(result.message()))
          .append(System.lineSeparator());
    }
    writeUtf8(output, csv.toString());
  }

  private static void writeUtf8(Path output, String content) throws IOException {
    Path parent = output.getParent();
    if (parent != null) {
      Files.createDirectories(parent);
    }
    Files.writeString(
        output,
        content,
        StandardCharsets.UTF_8,
        StandardOpenOption.CREATE,
        StandardOpenOption.TRUNCATE_EXISTING,
        StandardOpenOption.WRITE);
  }

  private static String format(double value) {
    return Double.isFinite(value) ? String.format(Locale.ROOT, "%.10f", value) : "";
  }

  private static String csv(String value) {
    String safe = value == null ? "" : value;
    return '"' + safe.replace("\"", "\"\"") + '"';
  }

  private static void validateDirectory(Path path, String label) throws IOException {
    if (path == null || !Files.isDirectory(path)) {
      throw new IOException("Invalid " + label + ": " + path);
    }
  }

  private static String rootMessage(Throwable throwable) {
    Throwable current = throwable;
    while (current.getCause() != null) {
      current = current.getCause();
    }
    String message = current.getMessage();
    return current.getClass().getSimpleName()
        + (message == null || message.isBlank() ? "" : ": " + message);
  }

  private static void printUsage() {
    System.err.println(
        "Usage: java ... org.tm_msaligner.TcoffeeBaselineMain "
            + "[precomputedRoot] [gpcrdbDataRoot] [consolidatedOutputCsv]");
  }

  private record Result(
      String dataset,
      String method,
      Path fasta,
      double objective1,
      double objective2,
      String status,
      String message) {

    static Result success(
        String dataset, String method, Path fasta, double objective1, double objective2) {
      return new Result(dataset, method, fasta, objective1, objective2, "ok", "");
    }

    static Result failure(String dataset, String method, Path fasta, String message) {
      return new Result(dataset, method, fasta, Double.NaN, Double.NaN, "failed", message);
    }
  }
}
