package org.tm_msaligner.evaluation;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class EvaluationMainTest {
  @TempDir
  Path tempDirectory;

  @Test
  void recursivelyDiscoversSeedAlignmentFilesInStableOrder() throws Exception {
    Path seed2 = Files.createDirectories(tempDirectory.resolve("seed_002"));
    Path seed1 = Files.createDirectories(tempDirectory.resolve("seed_001"));
    Files.writeString(seed2.resolve("MSA2.fasta"), ">a\nA\n");
    Files.writeString(seed1.resolve("MSA1.fasta"), ">a\nA\n");
    Files.writeString(seed1.resolve("notes.txt"), "ignored");

    List<Path> result = EvaluationMain.discoverAlignments(tempDirectory, "MSA*.fasta");

    assertEquals(List.of(
        seed1.resolve("MSA1.fasta"),
        seed2.resolve("MSA2.fasta")), result);
  }
}
