package org.tm_msaligner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.uma.jmetal.util.errorchecking.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

class TM_MSAStructAlignerGPCRdbTest {

  @TempDir
  Path temporaryDirectory;

  @Test
  void parsesAnExplicitRandomSeedAsTheNinthArgument() {
    String[] arguments = {
        "C:\\datasets\\GPCRdb",
        "classA_001",
        "2500",
        "100",
        "4",
        "100",
        "run_03",
        "C:\\results",
        "42"
    };

    TM_MSAStructAlignerGPCRdb.RunConfiguration configuration =
        TM_MSAStructAlignerGPCRdb.parseArguments(arguments);

    assertEquals("C:\\datasets\\GPCRdb", configuration.dataDirectory());
    assertEquals("classA_001", configuration.problemName());
    assertEquals(2500, configuration.maxEvaluations());
    assertEquals(100, configuration.populationSize());
    assertEquals(4, configuration.numberOfCores());
    assertEquals(100, configuration.observerFrequency());
    assertEquals("run_03", configuration.runId());
    assertEquals("C:\\results", configuration.outputRoot());
    assertEquals(42L, configuration.seed());
  }

  @Test
  void initializesTheGlobalJMetalGeneratorWithTheExplicitSeed() {
    TM_MSAStructAlignerGPCRdb.initializeRandomSeed(147L);
    assertEquals(147L, JMetalRandom.getInstance().getSeed());
  }

  @Test
  void rejectsTheOldEightArgumentInvocationWithoutASeed() {
    String[] arguments = {
        "C:\\datasets\\GPCRdb",
        "classA_001",
        "2500",
        "100",
        "4",
        "100",
        "run_01",
        "C:\\results"
    };

    JMetalException exception = assertThrows(
        JMetalException.class,
        () -> TM_MSAStructAlignerGPCRdb.parseArguments(arguments));

    assertTrue(exception.getMessage().contains("Expected 9"));
    assertTrue(exception.getMessage().contains("<seed>"));
  }

  @Test
  void writesSeedAndRuntimeMetadataUsingStableTextFormats() throws Exception {
    Files.createDirectories(temporaryDirectory);

    TM_MSAStructAlignerGPCRdb.writeSeedFile(temporaryDirectory, 126L);
    TM_MSAStructAlignerGPCRdb.writeRuntimeFile(
        temporaryDirectory,
        12345L,
        2500L);

    String seedText = Files.readString(
        temporaryDirectory.resolve("seed.txt"),
        StandardCharsets.UTF_8);

    String runtimeText = Files.readString(
        temporaryDirectory.resolve("runtime.txt"),
        StandardCharsets.UTF_8);

    assertEquals("126" + System.lineSeparator(), seedText);
    assertTrue(runtimeText.contains("runtime_ms\t12345"));
    assertTrue(runtimeText.contains("runtime_seconds\t12.345"));
    assertTrue(runtimeText.contains("evaluations\t2500"));
  }

  @Test
  void rejectsNonPositiveExperimentParameters() {
    String[] arguments = {
        "C:\\datasets\\GPCRdb",
        "classA_001",
        "0",
        "100",
        "4",
        "100",
        "run_01",
        "C:\\results",
        "0"
    };

    JMetalException exception = assertThrows(
        JMetalException.class,
        () -> TM_MSAStructAlignerGPCRdb.parseArguments(arguments));

    assertEquals("maxEvaluations must be greater than zero", exception.getMessage());
  }
}
