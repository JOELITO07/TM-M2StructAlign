package org.tm_msaligner.test;

import org.tm_msaligner.crossover.BioSPXMSACrossover;
import org.tm_msaligner.mutation.BioShiftClosedGapsMSAMutation;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.score.impl.LDDTStructuralScore;
import org.tm_msaligner.score.impl.StructSumOfPairsWithTopologyPredict;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AAArray;
import org.tm_msaligner.util.substitutionmatrix.impl.Blosum62;
import org.tm_msaligner.util.substitutionmatrix.impl.Phat;
import org.uma.jmetal.util.errorchecking.JMetalException;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Stress test for CROSSOVER (based on your mutation test example).
 *
 * - Loads GPCRdb problem + precomputed alignments
 * - Builds parents from precomputed alignments
 * - Applies crossover many times
 * - Validates offspring:
 *    1) isValid()
 *    2) alignment length consistent with (|seq| + #gaps) for each sequence
 *
 * Notes:
 * - Assumes BioSPXMSACrossover(probability, alpha) exists.
 *   If your operator still has only (probability), just adapt the constructor call.
 */
public class TM_MSAStructAlignerGPCRdb_CrossoverTest {

  public static void main(String[] args) throws JMetalException, IOException {

    // ===============================
    // 1) Parameters (same style as your example)
    // ===============================
    String dataDirectory = "D:\\Nube\\TM-MSA\\Datasets\\GPCRdb";
    String problemName = "classT2";
    String distanceDir = dataDirectory + "\\distances\\" + problemName;

    double probabilityCrossover = 0.8;
    double probabilityMutation = 0.2;

    // Bio-SPX scanning parameter (accept risky cut with prob = 1-alpha)
    double alpha = 0.8;

    int CROSSOVER_TESTS_PER_PAIR = 5000;

    var weightGapOpenTM = 8;
    var weightGapExtendTM = 3;
    var weightGapOpenNonTM = 3;
    var weightGapExtendNonTM = 1;

    String dataFile = dataDirectory + "\\sequences\\tmregions\\"
        + problemName + "_predicted_topologies.3line";

    // ===============================
    // 2) Load precomputed fasta alignments
    // ===============================
    List<String> preComputedFiles =
        getFastaFileNameListFromDir(dataDirectory + "\\precomputed\\" + problemName + "\\");

    if (preComputedFiles.size() < 2) {
      throw new JMetalException("Wrong number of Pre-computed Alignments, Minimum 2 files are required");
    }

    // ===============================
    // 3) Scores + Problem
    // ===============================
    List<StructuralScore> scoreList = new ArrayList<>();
    scoreList.add(new StructSumOfPairsWithTopologyPredict(
        new Phat(8),
        new Blosum62(),
        weightGapOpenTM,
        weightGapExtendTM,
        weightGapOpenNonTM,
        weightGapExtendNonTM
    ));
    scoreList.add(new LDDTStructuralScore(15.0f));

    MultiObjStructuralTMMSAProblem problem =
        new MultiObjStructuralTMMSAProblem(dataFile, scoreList, preComputedFiles, distanceDir, problemName);

    System.out.println("Problem loaded successfully!");
    System.out.println("Number of sequences: " + problem.numberOfVariables());

    // ===============================
    // 4) Operators (Mutation optional: only to diversify parents)
    // ===============================
    BioShiftClosedGapsMSAMutation mutationOperator = new BioShiftClosedGapsMSAMutation(probabilityMutation);

    // If your BioSPXMSACrossover DOESN'T have alpha yet, change to: new BioSPXMSACrossover(probabilityCrossover)
    BioSPXMSACrossover crossoverOperator = new BioSPXMSACrossover(probabilityCrossover, alpha);

    // ===============================
    // 5) Precomputed alignments loaded in memory (AAArray lists)
    // ===============================
    List<List<AAArray>> preAlignments = problem.listOfPrecomputedStringAlignments;
    System.out.println("Precomputed alignments loaded: " + preAlignments.size());

    long totalTests = 0;
    long totalFailures = 0;
    long totalExceptions = 0;

    // ===============================
    // 6) Test each pair of precomputed alignments (a vs b)
    //    - we test ALL pairs a<b (more coverage)
    // ===============================
    for (int a = 0; a < preAlignments.size(); a++) {
      for (int b = a + 1; b < preAlignments.size(); b++) {

        System.out.println("\n================================");
        System.out.println("Testing parents from alignments: A=" + a + " vs B=" + b);
        System.out.println("================================");

        StructuralTM_MSASolution parentA = new StructuralTM_MSASolution(preAlignments.get(a), problem);
        StructuralTM_MSASolution parentB = new StructuralTM_MSASolution(preAlignments.get(b), problem);

        if (!parentA.isValid()) {
          System.out.println("❌ ParentA INVALID. Skipping pair.");
          totalFailures++;
          continue;
        }
        if (!parentB.isValid()) {
          System.out.println("❌ ParentB INVALID. Skipping pair.");
          totalFailures++;
          continue;
        }

        // Important: crossover assumes SAME #sequences (variables size)
        if (parentA.variables().size() != parentB.variables().size()) {
          System.out.println("❌ Parents have different number of sequences. Skipping pair.");
          totalFailures++;
          continue;
        }

       
        // ===============================
        // 7) Stress test crossovers
        // ===============================
        for (int i = 0; i < CROSSOVER_TESTS_PER_PAIR; i++) {

          totalTests++;

          try {

            // Optionally diversify parents a bit each iteration (comment out if you want pure crossover)
            StructuralTM_MSASolution p1 = parentA.copy();
            StructuralTM_MSASolution p2 = parentB.copy();

            // light mutation to diversify parents (optional)
            // mutationOperator.execute(p1);
            // mutationOperator.execute(p2);

            // Build parents list for execute()
            List<StructuralTM_MSASolution> parents = new ArrayList<>(2);
            parents.add(p1);
            parents.add(p2);

            // Apply crossover
            List<StructuralTM_MSASolution> children = crossoverOperator.execute(parents);

            if (children == null || children.size() != 2) {
              System.out.println("❌ Crossover returned invalid children list at iter " + i);
              totalFailures++;
              continue;
            }

            // Validate children
            for (int c = 0; c < children.size(); c++) {
              StructuralTM_MSASolution child = children.get(c);

              if (!child.isValid()) {
                System.out.println("❌ Invalid child after crossover at iter " + i + " child=" + c);
                totalFailures++;
              }

              // Extra: alignment length consistency check per sequence
              int L = child.getAlignmentLength();
              for (int s = 0; s < child.variables().size(); s++) {
                int expected =
                    child.getOriginalSequences().get(s).getSize()
                        + child.getNumberOfGaps(s);

                if (L != expected) {
                  System.out.println("❌ Length mismatch after crossover | iter " + i
                      + " child=" + c + " seq=" + s
                      + " L=" + L + " expected=" + expected);
                  totalFailures++;
                  break;
                }
              }
            }

          } catch (Exception e) {
            totalExceptions++;
            totalFailures++;
            System.out.println("💥 Exception during crossover | A=" + a + " B=" + b + " iter " + i
                + " : " + e.getMessage());
            e.printStackTrace();
          }
        }

        System.out.println("Pair done A=" + a + " B=" + b
            + " | tests=" + CROSSOVER_TESTS_PER_PAIR
            + " | cumulative failures=" + totalFailures
            + " | exceptions=" + totalExceptions);
      }
    }

    // ===============================
    // 8) Final report
    // ===============================
    System.out.println("\n================================");
    System.out.println("Total crossover tests: " + totalTests);
    System.out.println("Total failures: " + totalFailures);
    System.out.println("Total exceptions: " + totalExceptions);
    System.out.println("================================");

    if (totalFailures == 0) {
      System.out.println("✅ All crossovers executed successfully.");
    } else {
      System.out.println("⚠ Some crossovers failed.");
    }
  }

  // Same helper you already have
  public static List<String> getFastaFileNameListFromDir(String dataDirectory) {
    List<String> preComputedFiles = new ArrayList<>();

    File File_Directory = new File(dataDirectory);
    if (!(File_Directory.exists() && File_Directory.isDirectory())) {
      System.out.println(String.format(dataDirectory + " does not exist"));
      return preComputedFiles;
    }
    FileFilter Demo_Filefilter = new FileFilter() {
      public boolean accept(File Demo_File) {
        return Demo_File.getName().endsWith(".fasta");
      }
    };

    File[] Text_Files = File_Directory.listFiles(Demo_Filefilter);
    if (Text_Files == null) return preComputedFiles;

    for (File Demo_File : Text_Files) {
      preComputedFiles.add(dataDirectory + Demo_File.getName());
    }

    return preComputedFiles;
  }
}