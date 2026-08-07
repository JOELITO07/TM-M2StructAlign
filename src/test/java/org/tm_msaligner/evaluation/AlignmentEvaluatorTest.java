package org.tm_msaligner.evaluation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class AlignmentEvaluatorTest {
  @TempDir
  Path tempDirectory;

  @Test
  void identicalAlignmentHasPerfectPairScoresDespiteOrderAndPdbExtension() throws Exception {
    Path topology = write("topology.3line", """
        >seq1
        ABCD
        OMMO
        >seq2
        ABCD
        OMMO
        """);
    Path reference = write("reference.fasta", """
        >seq1
        ABCD
        >seq2
        ABCD
        """);
    Path test = write("test.fasta", """
        >seq2.pdb
        ABCD
        >seq1.pdb:A
        ABCD
        """);

    AlignmentEvaluationResult result =
        new AlignmentEvaluator().evaluate(reference, test, topology);

    assertEquals(1.0, result.pairPrecision(), 1.0e-12);
    assertEquals(1.0, result.pairRecall(), 1.0e-12);
    assertEquals(1.0, result.pairF1(), 1.0e-12);
    assertEquals(1.0, result.tmPairF1(), 1.0e-12);
    assertEquals(0.0, result.tmGapRate(), 1.0e-12);
  }

  @Test
  void computesPartialPairF1AndTmGapRateFromUngappedResiduePositions() throws Exception {
    Path topology = write("topology.3line", """
        >seq1
        ABC
        OMM
        >seq2
        ABC
        OMM
        """);
    Path reference = write("reference.fasta", """
        >seq1
        AB-C
        >seq2
        A-BC
        """);
    Path test = write("test.fasta", """
        >seq1
        ABC-
        >seq2
        A-BC
        """);

    AlignmentEvaluationResult result =
        new AlignmentEvaluator().evaluate(reference, test, topology);

    assertEquals(2, result.referencePairCount());
    assertEquals(2, result.testPairCount());
    assertEquals(1, result.correctPairCount());
    assertEquals(0.5, result.pairPrecision(), 1.0e-12);
    assertEquals(0.5, result.pairRecall(), 1.0e-12);
    assertEquals(0.5, result.pairF1(), 1.0e-12);

    assertEquals(1, result.referenceTmPairCount());
    assertEquals(1, result.testTmPairCount());
    assertEquals(0, result.correctTmPairCount());
    assertEquals(0.0, result.tmPairF1(), 1.0e-12);

    assertEquals(2, result.tmGapCount());
    assertEquals(4, result.tmGapOpportunityCount());
    assertEquals(0.5, result.tmGapRate(), 1.0e-12);
  }

  @Test
  void rejectsASequenceMismatchInsteadOfReturningMisleadingScores() throws Exception {
    Path topology = write("topology.3line", """
        >seq1
        ABC
        MMM
        >seq2
        ABC
        MMM
        """);
    Path reference = write("reference.fasta", """
        >seq1
        ABC
        >seq2
        ABC
        """);
    Path test = write("test.fasta", """
        >seq1
        ABC
        >seq2
        ABD
        """);

    IOException exception = assertThrows(
        IOException.class,
        () -> new AlignmentEvaluator().evaluate(reference, test, topology));

    assertTrue(exception.getMessage().contains("Ungapped sequence mismatch"));
  }

  private Path write(String fileName, String content) throws IOException {
    Path path = tempDirectory.resolve(fileName);
    Files.writeString(path, content, StandardCharsets.UTF_8);
    return path;
  }
}
