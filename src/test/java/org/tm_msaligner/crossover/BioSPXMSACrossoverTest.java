package org.tm_msaligner.crossover;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.tm_msaligner.problem.StructuralTMMSAProblem;
import org.tm_msaligner.solution.StructuralTM_MSASolution;

class BioSPXMSACrossoverTest {

  @TempDir
  Path tempDirectory;

  @Test
  void detectsTheActualResidueBoundaryInsideAHelix() throws Exception {
    Fixture fixture = createFixture(
        "OMMMMO",
        List.of("AA---AAAA", "---AAAAAA"),
        List.of("A-A-AA-AA", "A--AAA-AA"));

    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1.0, 0.8);

    assertFalse(crossover.splitsTMHelixAfterResidue(fixture.parentA, 0, 1));
    assertTrue(crossover.splitsTMHelixAfterResidue(fixture.parentA, 0, 2));
    assertTrue(crossover.splitsTMHelixAfterResidue(fixture.parentA, 0, 4));
    assertFalse(crossover.splitsTMHelixAfterResidue(fixture.parentA, 0, 5));
  }

  @Test
  void rejectsAnAlignedColumnWhoseResidueBoundarySplitsAHelix() throws Exception {
    Fixture fixture = createFixture(
        "OMMMMO",
        List.of("AA---AAAA", "---AAAAAA"),
        List.of("A-A-AA-AA", "A--AAA-AA"));

    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1.0, 0.8);

    /*
     * At parentA column 3, sequence 0 contains a gap and sequence 1 contains
     * its first non-TM residue. A column-only risk check can therefore miss the
     * fact that sequence 0 is being cut after residue 2, inside an M-M helix
     * boundary. The residue-boundary check must reject it.
     */
    assertFalse(crossover.isCompatibleCut(fixture.parentA, fixture.parentB, 3));

    /* Both sequences are cut after residue 5: M-O, a valid helix boundary. */
    assertTrue(crossover.isCompatibleCut(fixture.parentA, fixture.parentB, 7));
    assertEquals(List.of(7),
        crossover.findCompatibleCutColumns(fixture.parentA, fixture.parentB));
  }

  @Test
  void validatesTheEquivalentBoundaryInBothParentDirections() throws Exception {
    Fixture fixture = createFixture(
        "OMMMMO",
        List.of("AA---AAAA", "---AAAAAA"),
        List.of("A-A-AA-AA", "A--AAA-AA"));

    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1.0, 0.8);

    assertEquals(List.of(7),
        crossover.findCompatibleCutColumns(fixture.parentA, fixture.parentB));
    assertEquals(List.of(7),
        crossover.findCompatibleCutColumns(fixture.parentB, fixture.parentA));
  }

  @Test
  void neverFallsBackToAnUnsafeRandomCut() throws Exception {
    Fixture fixture = createFixture(
        "MMMMMM",
        List.of("AAAAAA", "AAAAAA"),
        List.of("AAAAAA", "AAAAAA"));

    /*
     * alpha=0 made the previous implementation accept every risky cut. The new
     * strict implementation must copy the parents when no safe internal
     * boundary exists.
     */
    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1.0, 0.0);

    assertTrue(crossover.findCompatibleCutColumns(
        fixture.parentA, fixture.parentB).isEmpty());

    List<StructuralTM_MSASolution> children = crossover.execute(
        List.of(fixture.parentA, fixture.parentB));

    assertEquals(fixture.parentA.variables(), children.get(0).variables());
    assertEquals(fixture.parentB.variables(), children.get(1).variables());
    assertNotSame(fixture.parentA, children.get(0));
    assertNotSame(fixture.parentB, children.get(1));
  }

  @Test
  void generatedChildrenRemainValidAcrossRepeatedCrossovers() throws Exception {
    Fixture fixture = createFixture(
        "OMMMMO",
        List.of("AA---AAAA", "---AAAAAA"),
        List.of("A-A-AA-AA", "A--AAA-AA"));

    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1.0, 0.8);

    for (int iteration = 0; iteration < 250; iteration++) {
      List<StructuralTM_MSASolution> children = crossover.execute(
          List.of(fixture.parentA, fixture.parentB));

      assertEquals(2, children.size());
      for (StructuralTM_MSASolution child : children) {
        assertTrue(child.isValid(),
            "Invalid child generated at iteration " + iteration);

        int alignmentLength = child.getAlignmentLength();
        for (int sequence = 0; sequence < child.variables().size(); sequence++) {
          int expectedLength = child.getOriginalSequences().get(sequence).getSize()
              + child.getNumberOfGaps(sequence);
          assertEquals(expectedLength, alignmentLength,
              "Length mismatch at iteration " + iteration
                  + ", sequence " + sequence);
        }
      }
    }
  }

  private Fixture createFixture(
      String topology,
      List<String> firstAlignment,
      List<String> secondAlignment) throws Exception {

    if (firstAlignment.size() != 2 || secondAlignment.size() != 2) {
      throw new IllegalArgumentException("The fixture requires exactly two sequences");
    }

    Path datasetFile = tempDirectory.resolve("fixture_predicted_topologies.3line");
    Path distanceDirectory = Files.createDirectories(tempDirectory.resolve("distances"));
    Path firstFasta = tempDirectory.resolve("parentA.fasta");
    Path secondFasta = tempDirectory.resolve("parentB.fasta");

    String sequence = "AAAAAA";
    String dataset = ">seq1\n"
        + sequence + "\n"
        + topology + "\n"
        + ">seq2\n"
        + sequence + "\n"
        + topology + "\n";

    Files.writeString(datasetFile, dataset, StandardCharsets.UTF_8);
    writeFasta(firstFasta, firstAlignment);
    writeFasta(secondFasta, secondAlignment);
    writeStructuralFiles(distanceDirectory, "seq1", sequence.length());
    writeStructuralFiles(distanceDirectory, "seq2", sequence.length());

    StructuralTMMSAProblem problem = new StructuralTMMSAProblem(
        datasetFile.toString(),
        List.of(firstFasta.toString(), secondFasta.toString()),
        distanceDirectory.toString());

    StructuralTM_MSASolution parentA = new StructuralTM_MSASolution(
        problem.listOfPrecomputedStringAlignments.get(0),
        problem);

    StructuralTM_MSASolution parentB = new StructuralTM_MSASolution(
        problem.listOfPrecomputedStringAlignments.get(1),
        problem);

    assertTrue(parentA.isValid());
    assertTrue(parentB.isValid());

    return new Fixture(parentA, parentB);
  }

  private void writeFasta(Path path, List<String> alignment) throws IOException {
    String content = ">seq1\n"
        + alignment.get(0) + "\n"
        + ">seq2\n"
        + alignment.get(1) + "\n";
    Files.writeString(path, content, StandardCharsets.UTF_8);
  }

  private void writeStructuralFiles(
      Path distanceDirectory,
      String sequenceName,
      int length) throws IOException {

    List<String> matrixRows = new ArrayList<>();
    for (int row = 0; row < length; row++) {
      StringBuilder line = new StringBuilder();
      for (int column = 0; column < length; column++) {
        if (column > 0) {
          line.append(',');
        }
        line.append(Math.abs(row - column) * 3.0f);
      }
      matrixRows.add(line.toString());
    }

    Files.write(
        distanceDirectory.resolve(sequenceName + "_D.csv"),
        matrixRows,
        StandardCharsets.UTF_8);

    StringBuilder indices = new StringBuilder();
    for (int index = 1; index <= length; index++) {
      if (index > 1) {
        indices.append(System.lineSeparator());
      }
      indices.append(index);
    }

    Files.writeString(
        distanceDirectory.resolve(sequenceName + "_idx.csv"),
        indices.toString(),
        StandardCharsets.UTF_8);
  }

  private static final class Fixture {
    private final StructuralTM_MSASolution parentA;
    private final StructuralTM_MSASolution parentB;

    private Fixture(
        StructuralTM_MSASolution parentA,
        StructuralTM_MSASolution parentB) {
      this.parentA = parentA;
      this.parentB = parentB;
    }
  }
}
