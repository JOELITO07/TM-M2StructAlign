package org.tm_msaligner;


import org.tm_msaligner.algorithm.multiobjective.TM_M2AlignBuilder;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2Align;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2AlignBuilder;
import org.tm_msaligner.mutation.BioShiftClosedGapsMSAMutation;
import org.tm_msaligner.mutation.ShiftClosedGapsMSAMutation;
import org.tm_msaligner.problem.StandardTMMSAProblem;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.problem.impl.MultiObjTMMSAProblem;
import org.tm_msaligner.score.Score;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.score.impl.LDDTStructuralScore;
import org.tm_msaligner.score.impl.StructSumOfPairsWithTopologyPredict;
import org.tm_msaligner.score.impl.SumOfPairsWithTopologyPredict;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.util.AAArray;
import org.tm_msaligner.util.observer.FrontPlotTM_MSAObserver;
import org.tm_msaligner.util.observer.StructTM_MSAFitnessWriteFileObserver;
import org.tm_msaligner.util.observer.TM_MSAFitnessPlotObserver;
import org.tm_msaligner.util.observer.TM_MSAFitnessWriteFileObserver;
import org.tm_msaligner.util.substitutionmatrix.impl.Blosum62;
import org.tm_msaligner.util.substitutionmatrix.impl.Phat;
import org.tm_msaligner.util.visualization.MSAViewerHtmlMainPage;
import org.tm_msaligner.util.visualization.MSAViewerHtmlPage;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.errorchecking.JMetalException;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.observer.Observer;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TM_MSAStructAlignerGPCRdb extends AbstractAlgorithmRunner {

    public static void main(String[] args) throws JMetalException, IOException {

        //Parameters
        /*if (args.length != 7) {
            throw new JMetalException("Wrong number of arguments") ;
        }*/

        String dataDirectory = "D:\\Nube\\TM-MSA\\Datasets\\GPCRdb" ; // args[0]; // "data/gpcrdb/classA"
        String problemName = "classT2"; // args[1]; // "classA"
        Integer maxEvaluations = 250; //Integer.parseInt(args[2]);  //25000
        Integer populationSize = 10; //Integer.parseInt(args[3]); //100
        Integer numberOfCores = 1;   //Integer.parseInt(args[4]);   //1
        //0: Ninguno 1: FitnessWriteFileObserver, 2: FitnessPlotObserver y 3: FrontPlotTM_MSAObserve
        int observerType = 1; //Integer.parseInt(args[5]);
        int frequencyObserver = 20; //Integer.parseInt(args[6]);*/
        String distanceDir = dataDirectory + "\\distances\\" + problemName;

        //Algorithm  Parameters
        double probabilityCrossover=0.8;
        double probabilityMutation=0.2;
        var weightGapOpenTM = 8;
        var weightGapExtendTM = 3;
        var weightGapOpenNonTM = 3;
        var weightGapExtendNonTM = 1;

        String dataFile = dataDirectory + "\\sequences\\tmregions\\" + problemName +"_predicted_topologies.3line";
        String outputFolder = dataDirectory + "results" + System.currentTimeMillis() + "/";
        if (!new File(outputFolder).mkdirs()){
            throw new JMetalException("Error creating Output Directory " + outputFolder) ;
        }


        List<String> preComputedFiles = getFastaFileNameListFromDir(dataDirectory +"\\precomputed\\" + problemName + "\\");
        if(preComputedFiles.size()<2){
            throw new JMetalException("Wrong number of Pre-computed Alignments, Minimum 2 files are required") ;
        }
              
        List<StructuralScore> scoreList = new ArrayList<>();
        scoreList.add(new StructSumOfPairsWithTopologyPredict(
                new Phat(8),
                new Blosum62(),
                weightGapOpenTM,
                weightGapExtendTM,
                weightGapOpenNonTM,
                weightGapExtendNonTM
            )
        );
        scoreList.add(new LDDTStructuralScore(  15.0f)     );


        MultiObjStructuralTMMSAProblem problem = new MultiObjStructuralTMMSAProblem(dataFile, scoreList,
                                                                 preComputedFiles,  distanceDir, problemName);

        System.out.println("Problem loaded successfully!");
        System.out.println("Number of sequences: " + problem.numberOfVariables());
        
        BioShiftClosedGapsMSAMutation mutationOperator =  new BioShiftClosedGapsMSAMutation(probabilityMutation); 

        /*try {
        
               

                List<List<AAArray>> preAlignments = problem.listOfPrecomputedStringAlignments;
                System.out.println("Precomputed alignments loaded: " + preAlignments.size());
                int totalTests = 0;
                int totalFailures = 0;

                // ===============================
                // 4️⃣  Testear cada alineamiento
                // ===============================
                for (int a = 0; a < preAlignments.size(); a++) {

                    System.out.println("\nTesting alignment: " + a);

                    StructuralTM_MSASolution solution =  new StructuralTM_MSASolution(preAlignments.get(a), problem);

                    if (!solution.isValid()) {
                        System.out.println("Initial solution INVALID!");
                        continue;
                    }

                        // ===============================
                        // 5️⃣  Aplicar muchas mutaciones
                        // ===============================
                        for (int i = 0; i < 5000; i++) {

                        totalTests++;

                        try {

                            StructuralTM_MSASolution mutated = solution.copy();

                            mutationOperator.execute(mutated);

                            if (!mutated.isValid()) {
                                System.out.println("❌ Invalid solution after mutation at iteration " + i);
                                totalFailures++;
                            }

                            // Validación extra: longitud consistente
                            int L = mutated.getAlignmentLength();
                            for (int s = 0; s < mutated.variables().size(); s++) {
                            int expected =
                                    mutated.getOriginalSequences().get(s).getSize()
                                            + mutated.getNumberOfGaps(s);

                            if (L != expected) {
                                System.out.println("❌ Length mismatch after mutation");
                                totalFailures++;
                            }
                            }

                        } catch (Exception e) {
                            System.out.println("💥 Exception during mutation: " + e.getMessage());
                            e.printStackTrace();
                            totalFailures++;
                        }
                        }
                    }

            // ===============================
            // 6️⃣  Reporte final
            // ===============================
            System.out.println("\n================================");
            System.out.println("Total mutation tests: " + totalTests);
            System.out.println("Total failures: " + totalFailures);
            System.out.println("================================");

            if (totalFailures == 0) {
                System.out.println("✅ All mutations executed successfully.");
            } else {
                System.out.println("⚠ Some mutations failed.");
            }

        } catch (Exception e) {
          e.printStackTrace();
        }*/

   
        int offspringPopulationSize = populationSize;
        StructuralTM_M2Align tm_m2align = new StructuralTM_M2AlignBuilder(problem,
                            maxEvaluations,
                            populationSize,
                            offspringPopulationSize,
                            probabilityCrossover,
                            mutationOperator,
                            numberOfCores)
                            .build();


        if(observerType>=1 && observerType<=3){

            if(frequencyObserver> maxEvaluations){
                throw new JMetalException("The frequency of the Observer can`t be greater than Maximun number of Evaluations") ;
            }

            Observer chartObserver;
            chartObserver = new StructTM_MSAFitnessWriteFileObserver(outputFolder + "BestScores.tsv", frequencyObserver);
           
           /* if(observerType==1) {
                chartObserver = new TM_MSAFitnessWriteFileObserver(outputFolder + "BestScores.tsv", 100);
            } else if (observerType==2) {
               // chartObserver = new TM_MSAFitnessPlotObserver("TM-M2Aligner solving Instance " + problemName ,
               //         "Evaluations", scoreList.get(0).getName(), scoreList.get(0).getName(), frequencyObserver, 0);
            }else{
               // chartObserver = new FrontPlotTM_MSAObserver<TM_MSASolution>("", "SumOfPairsWithTopologyPredict",
               //         "AlignedSegment", problemName, frequencyObserver);
            }*/

            tm_m2align.observable().register(chartObserver);
        }


        tm_m2align.run();
        List<StructuralTM_MSASolution> population = tm_m2align.result();

        for (StructuralTM_MSASolution solution : population) {
            for (int i = 0; i < problem.numberOfObjectives(); i++) {
                solution.objectives()[i] *= (scoreList.get(i).isAMinimizationScore()?1.0:-1.0);
            }
        }

        JMetalLogger.logger.info("Total execution time : " + tm_m2align.totalComputingTime() + "ms");
        JMetalLogger.logger.info("Number of evaluations: " + tm_m2align.numberOfEvaluations()) ;

        DefaultFileOutputContext funFile = new DefaultFileOutputContext(outputFolder + "FUN.tsv");
        funFile.setSeparator("\t");

        SolutionListOutput slo = new SolutionListOutput(population);
        slo.printObjectivesToFile(funFile, population);

        printMSASolutionsToFile(population, outputFolder);
        




    }

    public static List<String> getFastaFileNameListFromDir(String dataDirectory){
        List<String> preComputedFiles = new ArrayList<>();

        File File_Directory = new File(dataDirectory);
        if (!(File_Directory.exists() && File_Directory.isDirectory())) {
            System.out.println(String.format(dataDirectory + " does not exist"));
            return preComputedFiles;
        }
        FileFilter Demo_Filefilter = new FileFilter() {
            public boolean accept(File Demo_File) {
                if (Demo_File.getName().endsWith(".fasta")) return true;
                return false;
            }
        };

        File[] Text_Files = File_Directory.listFiles(Demo_Filefilter);
        for (File Demo_File: Text_Files)
            preComputedFiles.add(dataDirectory + Demo_File.getName());


        return preComputedFiles;
    }

    public static void printMSASolutionsToFile(List<StructuralTM_MSASolution> solutionList, String PathOut) {
        for (int i = 0; i < solutionList.size(); i++)
            solutionList.get(i).printSolutionToFasta(PathOut + "MSASol" + i + ".fasta");
    }
}
