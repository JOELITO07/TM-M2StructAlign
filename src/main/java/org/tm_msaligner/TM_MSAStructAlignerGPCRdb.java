package org.tm_msaligner;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2Align;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2AlignBuilder;
import org.tm_msaligner.crossover.BioSPXMSACrossover;
import org.tm_msaligner.mutation.BioShiftClosedGapsMSAMutation;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.score.impl.LDDTStructuralScore;
import org.tm_msaligner.score.impl.StructSumOfPairsWithTopologyPredict;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.observer.StructTM_MSAFitnessWriteFileObserver;
import org.tm_msaligner.util.substitutionmatrix.impl.Blosum62;
import org.tm_msaligner.util.substitutionmatrix.impl.Phat;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.errorchecking.JMetalException;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.observer.Observer;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

/** Runner for the structural multiobjective GPCRdb experiments. */
public class TM_MSAStructAlignerGPCRdb extends AbstractAlgorithmRunner {

  private static final String USAGE =
      "Usage: java -jar TM-MSAStructAligner.jar "
          + "<dataDirectory> <problemName> <maxEvaluations> <populationSize> "
          + "<numberOfCores> <observerFrequency> <runId> <outputRoot> <seed>";

  public static void main(String[] args) throws JMetalException, IOException {
    RunConfiguration configuration = parseArguments(args);

    /*
     * Set the global jMetal seed before creating operators, builders, or the
     * initial population. Every stochastic component using JMetalRandom will
     * therefore consume the same pseudo-random sequence when the same seed and
     * configuration are used.
     */
    initializeRandomSeed(configuration.seed());

    Path distanceDirectory = Paths.get(
        configuration.dataDirectory(),
        "distances",
        configuration.problemName());

    Path dataFilePath = Paths.get(
        configuration.dataDirectory(),
        "sequences",
        "tmregions",
        configuration.problemName() + "_predicted_topologies.3line");

    Path outputFolder = Paths.get(
        configuration.outputRoot(),
        "ejecuciones",
        configuration.problemName(),
        configuration.runId());

    Path precomputedFolder = Paths.get(
        configuration.dataDirectory(),
        "precomputed",
        configuration.problemName());

    if (Files.exists(outputFolder)) {
      throw new JMetalException(
          "Output directory already exists and will not be overwritten: "
              + outputFolder);
    }
    Files.createDirectories(outputFolder);

    /* The seed is persisted before the algorithm starts. */
    writeSeedFile(outputFolder, configuration.seed());

    List<String> precomputedFiles = getFastaFileNameListFromDir(
        precomputedFolder.toString());
    if (precomputedFiles.size() < 2) {
      throw new JMetalException(
          "Wrong number of Pre-computed Alignments, Minimum 2 files are required");
    }

    double probabilityCrossover = 0.8;
    double probabilityMutation = 0.2;
    double alpha = 0.8;
    int weightGapOpenTM = 8;
    int weightGapExtendTM = 3;
    int weightGapOpenNonTM = 3;
    int weightGapExtendNonTM = 1;

    List<StructuralScore> scoreList = new ArrayList<>();
    scoreList.add(
        new StructSumOfPairsWithTopologyPredict(
            new Phat(8),
            new Blosum62(),
            weightGapOpenTM,
            weightGapExtendTM,
            weightGapOpenNonTM,
            weightGapExtendNonTM));
    scoreList.add(new LDDTStructuralScore(15.0f));

    MultiObjStructuralTMMSAProblem problem = new MultiObjStructuralTMMSAProblem(
        dataFilePath.toString(),
        scoreList,
        precomputedFiles,
        distanceDirectory.toString(),
        configuration.problemName());

    System.out.println("Problem loaded successfully!");
    System.out.println("Number of sequences: " + problem.numberOfVariables());
    System.out.println("Random seed: " + configuration.seed());
    System.out.println("Run ID: " + configuration.runId());

    BioShiftClosedGapsMSAMutation mutationOperator =
        new BioShiftClosedGapsMSAMutation(probabilityMutation);
    BioSPXMSACrossover crossoverOperator =
        new BioSPXMSACrossover(probabilityCrossover, alpha);

    int offspringPopulationSize = configuration.populationSize();
    StructuralTM_M2Align tmM2Align = new StructuralTM_M2AlignBuilder(
        problem,
        configuration.maxEvaluations(),
        configuration.populationSize(),
        offspringPopulationSize,
        crossoverOperator,
        mutationOperator,
        configuration.numberOfCores())
        .build();

    Observer chartObserver = new StructTM_MSAFitnessWriteFileObserver(
        outputFolder.resolve("BestScores.tsv").toString(),
        configuration.observerFrequency());
    tmM2Align.observable().register(chartObserver);

    tmM2Align.run();
    List<StructuralTM_MSASolution> population = tmM2Align.result();

    for (StructuralTM_MSASolution solution : population) {
      for (int objective = 0; objective < problem.numberOfObjectives(); objective++) {
        solution.objectives()[objective] *=
            scoreList.get(objective).isAMinimizationScore() ? 1.0 : -1.0;
      }
    }

    long runtimeMilliseconds = tmM2Align.totalComputingTime();
    JMetalLogger.logger.info(
        "Total execution time: " + runtimeMilliseconds + " ms");
    JMetalLogger.logger.info(
        "Number of evaluations: " + tmM2Align.numberOfEvaluations());
    JMetalLogger.logger.info("Random seed: " + configuration.seed());

    DefaultFileOutputContext funFile = new DefaultFileOutputContext(
        outputFolder.resolve("FUN.tsv").toString());
    funFile.setSeparator("\t");

    SolutionListOutput solutionOutput = new SolutionListOutput(population);
    solutionOutput.printObjectivesToFile(funFile, population);

    printMSASolutionsToFile(population, outputFolder.toString());
    writeRuntimeFile(
        outputFolder,
        runtimeMilliseconds,
        tmM2Align.numberOfEvaluations());
  }

  static RunConfiguration parseArguments(String[] args) {
    if (args == null || args.length != 9) {
      throw new JMetalException(
          "Wrong number of arguments. Expected 9 but received "
              + (args == null ? 0 : args.length)
              + ".\n"
              + USAGE);
    }

    try {
      return new RunConfiguration(
          args[0],
          args[1],
          Integer.parseInt(args[2]),
          Integer.parseInt(args[3]),
          Integer.parseInt(args[4]),
          Integer.parseInt(args[5]),
          args[6],
          args[7],
          Long.parseLong(args[8]));
    } catch (NumberFormatException exception) {
      throw new JMetalException(
          "Invalid numeric argument: " + exception.getMessage() + ".\n" + USAGE);
    }
  }

  static void initializeRandomSeed(long seed) {
    JMetalRandom.getInstance().setSeed(seed);
  }

  static void writeSeedFile(Path outputFolder, long seed) throws IOException {
    Files.writeString(
        outputFolder.resolve("seed.txt"),
        Long.toString(seed) + System.lineSeparator(),
        StandardCharsets.UTF_8);
  }

  static void writeRuntimeFile(
      Path outputFolder,
      long runtimeMilliseconds,
      long evaluations) throws IOException {

    String content = "runtime_ms\t" + runtimeMilliseconds + System.lineSeparator()
        + "runtime_seconds\t"
        + String.format(Locale.ROOT, "%.3f", runtimeMilliseconds / 1000.0)
        + System.lineSeparator()
        + "evaluations\t"
        + evaluations
        + System.lineSeparator();

    Files.writeString(
        outputFolder.resolve("runtime.txt"),
        content,
        StandardCharsets.UTF_8);
  }

  public static List<String> getFastaFileNameListFromDir(String dataDirectory) {
    List<String> precomputedFiles = new ArrayList<>();

    File directory = new File(dataDirectory);
    if (!(directory.exists() && directory.isDirectory())) {
      System.out.println(dataDirectory + " does not exist");
      return precomputedFiles;
    }

    FileFilter fastaFilter = file -> file.getName().toLowerCase(Locale.ROOT).endsWith(".fasta");
    File[] fastaFiles = directory.listFiles(fastaFilter);

    if (fastaFiles == null) {
      return precomputedFiles;
    }

    java.util.Arrays.sort(fastaFiles, java.util.Comparator.comparing(File::getName));
    for (File fastaFile : fastaFiles) {
      precomputedFiles.add(Paths.get(dataDirectory, fastaFile.getName()).toString());
    }

    return precomputedFiles;
  }

  public static void printMSASolutionsToFile(
      List<StructuralTM_MSASolution> solutionList,
      String outputPath) {

    Path outputDirectory = Paths.get(outputPath);
    for (int index = 0; index < solutionList.size(); index++) {
      Path outputFile = outputDirectory.resolve("MSASol" + index + ".fasta");
      solutionList.get(index).printSolutionToFasta(outputFile.toString());
    }
  }

  static record RunConfiguration(
      String dataDirectory,
      String problemName,
      int maxEvaluations,
      int populationSize,
      int numberOfCores,
      int observerFrequency,
      String runId,
      String outputRoot,
      long seed) {

    RunConfiguration {
      if (dataDirectory == null || dataDirectory.isBlank()) {
        throw new JMetalException("dataDirectory must not be blank");
      }
      if (problemName == null || problemName.isBlank()) {
        throw new JMetalException("problemName must not be blank");
      }
      if (maxEvaluations <= 0) {
        throw new JMetalException("maxEvaluations must be greater than zero");
      }
      if (populationSize <= 0) {
        throw new JMetalException("populationSize must be greater than zero");
      }
      if (numberOfCores <= 0) {
        throw new JMetalException("numberOfCores must be greater than zero");
      }
      if (observerFrequency <= 0) {
        throw new JMetalException("observerFrequency must be greater than zero");
      }
      if (runId == null || runId.isBlank()) {
        throw new JMetalException("runId must not be blank");
      }
      if (outputRoot == null || outputRoot.isBlank()) {
        throw new JMetalException("outputRoot must not be blank");
      }
    }
  }
}
