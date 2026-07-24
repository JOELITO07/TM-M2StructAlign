package org.tm_msaligner;

import java.nio.file.Path;
import java.nio.file.Paths;
import org.tm_msaligner.algorithm.multiobjective.TM_M2AlignBuilder;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2Align;
import org.tm_msaligner.algorithm.structural_multiobjective.StructuralTM_M2AlignBuilder;
import org.tm_msaligner.crossover.BioSPXMSACrossover;
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
        if (args.length != 8) {
            throw new JMetalException("Wrong number of arguments") ;
        }

        String dataDirectory = args[0]; // "D:\\Nube\\TM-MSA\\Datasets\\GPCRdb" ; //  // "data/gpcrdb/classA"
        String problemName = args[1]; // "classF"; //  // "classA"
        Integer maxEvaluations = Integer.parseInt(args[2]); // 2500; //Integer.parseInt(args[2]);  //25000
        Integer populationSize = Integer.parseInt(args[3]); //100
        Integer numberOfCores = Integer.parseInt(args[4]);   //1
        int frequencyObserver = Integer.parseInt(args[5]); //Integer.parseInt(args[6]);*/
        String numberTest = args[6]; //Integer.parseInt(args[6]);*/
        String outputFolderPath = args[7];
        
        Path disPath = Paths.get(dataDirectory,"distances", problemName);
        String distanceDir = disPath.toString();

        //Algorithm  Parameters
        double probabilityCrossover=0.8;
        double probabilityMutation=0.2;
        double alpha = 0.8; 
        var weightGapOpenTM = 8;
        var weightGapExtendTM = 3;
        var weightGapOpenNonTM = 3;
        var weightGapExtendNonTM = 1;

        Path dataFilePath = Paths.get(dataDirectory,"sequences","tmregions", problemName + "_predicted_topologies.3line");
        Path outputFolder = Paths.get(outputFolderPath, "ejecuciones", problemName, numberTest );
        Path precomputedFolder = Paths.get(dataDirectory, "precomputed", problemName);


        List<String> preComputedFiles = getFastaFileNameListFromDir(precomputedFolder.toString());
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


        MultiObjStructuralTMMSAProblem problem = new MultiObjStructuralTMMSAProblem(dataFilePath.toString(), scoreList,
                                                                 preComputedFiles,  distanceDir, problemName);

        System.out.println("Problem loaded successfully!");
        System.out.println("Number of sequences: " + problem.numberOfVariables());
        
        BioShiftClosedGapsMSAMutation mutationOperator =  new BioShiftClosedGapsMSAMutation(probabilityMutation); 
        BioSPXMSACrossover crossoverOperator = new BioSPXMSACrossover(probabilityCrossover, alpha);
    
   
        int offspringPopulationSize = populationSize;
        StructuralTM_M2Align tm_m2align = new StructuralTM_M2AlignBuilder(problem,
                            maxEvaluations,
                            populationSize,
                            offspringPopulationSize,
                            crossoverOperator,
                            mutationOperator,
                            numberOfCores)
                            .build();

        
        if (!new File(outputFolder.toString()).mkdirs()){
            throw new JMetalException("Error creating Output Directory " + outputFolder.toString()) ;
        }

        Observer chartObserver;
        chartObserver = new StructTM_MSAFitnessWriteFileObserver(outputFolder.resolve("BestScores.tsv").toString(), frequencyObserver);
        tm_m2align.observable().register(chartObserver);

       


       tm_m2align.run();
        List<StructuralTM_MSASolution> population = tm_m2align.result();

        for (StructuralTM_MSASolution solution : population) {
            for (int i = 0; i < problem.numberOfObjectives(); i++) {
                solution.objectives()[i] *= (scoreList.get(i).isAMinimizationScore()?1.0:-1.0);
            }
        }

        JMetalLogger.logger.info("Total execution time : " + tm_m2align.totalComputingTime() + "ms");
        JMetalLogger.logger.info("Number of evaluations: " + tm_m2align.numberOfEvaluations()) ;

        DefaultFileOutputContext funFile = new DefaultFileOutputContext(outputFolder.resolve("FUN.tsv").toString());
        funFile.setSeparator("\t");

        SolutionListOutput slo = new SolutionListOutput(population);
        slo.printObjectivesToFile(funFile, population);

        printMSASolutionsToFile(population, outputFolder.toString());
        




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
            preComputedFiles.add(Paths.get(dataDirectory, Demo_File.getName()).toString());
        


        return preComputedFiles;
    }

    public static void printMSASolutionsToFile(List<StructuralTM_MSASolution> solutionList, String PathOut) {

        Path outDir = Paths.get(PathOut);
        for (int i = 0; i < solutionList.size(); i++){
            Path outFile = outDir.resolve("MSASol" + i + ".fasta");
            solutionList.get(i).printSolutionToFasta(outFile.toString());
        }
    }
}
