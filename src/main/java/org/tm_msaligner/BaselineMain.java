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

/** Evaluates every precomputed FASTA with the two project objectives. */
public final class BaselineMain {

  private static final Path DEFAULT_PRECOMPUTED_ROOT = Path.of(
      "C:\\GDrive2026\\TM-MSA\\Datasets\\GPCRdb\\precomputed");

  private BaselineMain() {
  }

  public static void main(String[] args) throws IOException {
    if (args.length > 3) {
      System.err.println("Usage: BaselineMain [precomputedRoot] [gpcrdbRoot] [outputCsv]");
      System.exit(2);
    }

    Path precomputedRoot = args.length >= 1
        ? Path.of(args[0]).toAbsolutePath().normalize()
        : DEFAULT_PRECOMPUTED_ROOT.toAbsolutePath().normalize();
    Path gpcrdbRoot = args.length >= 2
        ? Path.of(args[1]).toAbsolutePath().normalize()
        : precomputedRoot.getParent();
    Path output = args.length >= 3
        ? Path.of(args[2]).toAbsolutePath().normalize()
        : precomputedRoot.resolve("baselines.csv");

    requireDirectory(precomputedRoot, "precomputed root");
    requireDirectory(gpcrdbRoot, "GPCRdb root");

    List<Path> datasets = directories(precomputedRoot);
    List<Result> allResults = new ArrayList<>();
    int successful = 0;
    int attempted = 0;

    System.out.println("Precomputed: " + precomputedRoot);
    System.out.println("GPCRdb root: " + gpcrdbRoot);
    System.out.println("Datasets: " + datasets.size());

    for (Path datasetDirectory : datasets) {
      String dataset = datasetDirectory.getFileName().toString();
      List<Path> fastaFiles = fastaFiles(datasetDirectory);
      List<Result> datasetResults = new ArrayList<>();

      System.out.println();
      System.out.println("Dataset " + dataset + " (" + fastaFiles.size() + " FASTA)");

      for (Path fasta : fastaFiles) {
        attempted++;
        String method = methodName(fasta, dataset);
        try {
          Result result = evaluate(gpcrdbRoot, dataset, method, fasta);
          datasetResults.add(result);
          allResults.add(result);
          successful++;
          System.out.printf(Locale.ROOT,
              "  [OK] %-20s objective_1=%.10f objective_2=%.10f%n",
              method, result.objective1(), result.objective2());
        } catch (Exception exception) {
          Result result = Result.failure(dataset, method, fasta, rootMessage(exception));
          datasetResults.add(result);
          allResults.add(result);
          System.err.println("  [FAIL] " + method + ": " + result.message());
        }
      }

      // This is the file consumed by pareto_analysis.py --baselines.
      writeDatasetCsv(datasetDirectory.resolve("baselines.csv"), datasetResults);
    }

    writeConsolidatedCsv(output, allResults);
    System.out.println();
    System.out.println("Successful alignments: " + successful + "/" + attempted);
    System.out.println("Consolidated output: " + output);

    if (successful == 0) {
      System.exit(1);
    }
  }

  private static Result evaluate(Path gpcrdbRoot, String dataset, String method, Path fasta)
      throws IOException {
    Path topology = gpcrdbRoot.resolve("sequences").resolve("tmregions")
        .resolve(dataset + "_predicted_topologies.3line");
    Path distances = gpcrdbRoot.resolve("distances").resolve(dataset);

    if (!Files.isRegularFile(topology)) {
      throw new IOException("Missing topology file: " + topology);
    }
    requireDirectory(distances, "distance directory for " + dataset);

    List<StructuralScore> scores = scores();

    /* The problem constructor requires two initial alignments. Repeating the
       target only satisfies construction; no evolutionary run is performed. */
    MultiObjStructuralTMMSAProblem problem = new MultiObjStructuralTMMSAProblem(
        topology.toString(), scores, List.of(fasta.toString(), fasta.toString()),
        distances.toString(), dataset);

    List<AAArray> aligned = problem.readDataFromFastaFile(fasta.toString());
    StructuralTM_MSASolution solution = new StructuralTM_MSASolution(aligned, problem);
    if (!solution.isValid()) {
      throw new IllegalArgumentException("Invalid MSA");
    }

    // Same preprocessing performed by MultiObjStructuralTMMSAProblem.evaluate().
    solution.removeGapColumns();
    AA[][] decoded = solution.decodeToMatrix();

    // Raw benefit values, identical to those written to FUN.tsv by the runner.
    double objective1 = scores.get(0).compute(solution, decoded);
    double objective2 = scores.get(1).compute(solution, decoded);
    if (!Double.isFinite(objective1) || !Double.isFinite(objective2)) {
      throw new IllegalStateException("Non-finite objectives");
    }
    return Result.success(dataset, method, fasta, objective1, objective2);
  }

  private static List<StructuralScore> scores() {
    return List.of(
        new StructSumOfPairsWithTopologyPredict(
            new Phat(8), new Blosum62(), 8.0, 3.0, 3.0, 1.0),
        new LDDTStructuralScore(15.0f));
  }

  private static List<Path> directories(Path root) throws IOException {
    List<Path> result = new ArrayList<>();
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(root)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          result.add(path);
        }
      }
    }
    result.sort(Comparator.comparing(path -> path.getFileName().toString(),
        String.CASE_INSENSITIVE_ORDER));
    return result;
  }

  private static List<Path> fastaFiles(Path directory) throws IOException {
    List<Path> result = new ArrayList<>();
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(directory)) {
      for (Path path : stream) {
        if (Files.isRegularFile(path)
            && path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".fasta")) {
          result.add(path);
        }
      }
    }
    result.sort(Comparator.comparing(path -> path.getFileName().toString(),
        String.CASE_INSENSITIVE_ORDER));
    return result;
  }

  private static String methodName(Path fasta, String dataset) {
    String stem = fasta.getFileName().toString();
    stem = stem.substring(0, stem.length() - ".fasta".length());
    String token = "_" + dataset;
    int position = stem.toLowerCase(Locale.ROOT).indexOf(token.toLowerCase(Locale.ROOT));
    if (position < 0) {
      return stem;
    }
    String prefix = stem.substring(0, position);
    String suffix = stem.substring(position + token.length());
    return prefix + suffix; // e.g. usalign_classA_001_clean -> usalign_clean
  }

  private static void writeDatasetCsv(Path output, List<Result> results) throws IOException {
    StringBuilder text = new StringBuilder("method,objective_1,objective_2")
        .append(System.lineSeparator());
    for (Result result : results) {
      if (result.ok()) {
        text.append(csv(result.method())).append(',')
            .append(number(result.objective1())).append(',')
            .append(number(result.objective2())).append(System.lineSeparator());
      }
    }
    write(output, text.toString());
  }

  private static void writeConsolidatedCsv(Path output, List<Result> results) throws IOException {
    StringBuilder text = new StringBuilder(
        "dataset,method,objective_1,objective_2,fasta,status,message")
        .append(System.lineSeparator());
    for (Result result : results) {
      text.append(csv(result.dataset())).append(',')
          .append(csv(result.method())).append(',')
          .append(number(result.objective1())).append(',')
          .append(number(result.objective2())).append(',')
          .append(csv(result.fasta().toString())).append(',')
          .append(result.ok() ? "ok" : "failed").append(',')
          .append(csv(result.message())).append(System.lineSeparator());
    }
    write(output, text.toString());
  }

  private static void write(Path output, String text) throws IOException {
    if (output.getParent() != null) {
      Files.createDirectories(output.getParent());
    }
    Files.writeString(output, text, StandardCharsets.UTF_8,
        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING,
        StandardOpenOption.WRITE);
  }

  private static String number(double value) {
    return Double.isFinite(value) ? String.format(Locale.ROOT, "%.10f", value) : "";
  }

  private static String csv(String value) {
    String safe = value == null ? "" : value;
    return "\"" + safe.replace("\"", "\"\"") + "\"";
  }

  private static void requireDirectory(Path path, String label) throws IOException {
    if (path == null || !Files.isDirectory(path)) {
      throw new IOException("Invalid " + label + ": " + path);
    }
  }

  private static String rootMessage(Throwable throwable) {
    Throwable current = throwable;
    while (current.getCause() != null) {
      current = current.getCause();
    }
    return current.getClass().getSimpleName() + ": "
        + (current.getMessage() == null ? "" : current.getMessage());
  }

  private record Result(String dataset, String method, Path fasta,
                        double objective1, double objective2,
                        boolean ok, String message) {
    static Result success(String dataset, String method, Path fasta,
                          double objective1, double objective2) {
      return new Result(dataset, method, fasta, objective1, objective2, true, "");
    }

    static Result failure(String dataset, String method, Path fasta, String message) {
      return new Result(dataset, method, fasta, Double.NaN, Double.NaN, false, message);
    }
  }
}
