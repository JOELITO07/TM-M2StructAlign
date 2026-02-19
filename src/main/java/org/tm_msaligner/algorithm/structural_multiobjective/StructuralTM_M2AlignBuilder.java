package org.tm_msaligner.algorithm.structural_multiobjective;

import java.util.Arrays;
import java.util.Comparator;

import org.tm_msaligner.algorithm.multiobjective.TM_M2Align;
import org.tm_msaligner.crossover.BioSPXMSACrossover;
import org.tm_msaligner.crossover.SPXMSACrossover;
import org.tm_msaligner.mutation.ShiftClosedGapsMSAMutation;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.solutionscreation.PreComputedMSAsSolutionsCreation;
import org.tm_msaligner.solutionscreation.PreComputedStructMSAsSolutionsCreation;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.MultiThreadedEvaluation;
import org.uma.jmetal.component.catalogue.common.evaluation.impl.SequentialEvaluation;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.component.catalogue.common.termination.impl.TerminationByEvaluations;
import org.uma.jmetal.component.catalogue.ea.replacement.Replacement;
import org.uma.jmetal.component.catalogue.ea.replacement.impl.RankingAndDensityEstimatorReplacement;
import org.uma.jmetal.component.catalogue.ea.selection.Selection;
import org.uma.jmetal.component.catalogue.ea.selection.impl.NaryTournamentSelection;
import org.uma.jmetal.component.catalogue.ea.variation.Variation;
import org.uma.jmetal.component.catalogue.ea.variation.impl.CrossoverAndMutationVariation;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.operator.mutation.MutationOperator;
import org.uma.jmetal.util.comparator.MultiComparator;
import org.uma.jmetal.util.densityestimator.DensityEstimator;
import org.uma.jmetal.util.densityestimator.impl.CrowdingDistanceDensityEstimator;
import org.uma.jmetal.util.ranking.Ranking;
import org.uma.jmetal.util.ranking.impl.FastNonDominatedSortRanking;
import org.tm_msaligner.problem.StandardTMMSAProblem;
import org.tm_msaligner.problem.StructuralTMMSAProblem;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;

public class StructuralTM_M2AlignBuilder {

  Ranking<StructuralTM_MSASolution> ranking;
  private Evaluation<StructuralTM_MSASolution> evaluation;
  private PreComputedStructMSAsSolutionsCreation createInitialPopulation;
  private Termination termination;
  private Selection<StructuralTM_MSASolution> selection;
  private Variation<StructuralTM_MSASolution> variation;
  private Replacement<StructuralTM_MSASolution> replacement;
  private DensityEstimator<StructuralTM_MSASolution> densityEstimator;
  private CrossoverOperator<StructuralTM_MSASolution> crossover;
  private MutationOperator<StructuralTM_MSASolution> mutation;


  public StructuralTM_M2AlignBuilder(StructuralTMMSAProblem problem, int maxEvaluations,
      int populationSize, int offspringPopulationSize,
      double probabilityCrossover, MutationOperator<StructuralTM_MSASolution> mutationOperator,
      int numCores) {

    crossover = new BioSPXMSACrossover(probabilityCrossover);
    mutation = mutationOperator;
    variation = new CrossoverAndMutationVariation<>(
        offspringPopulationSize, crossover, mutation);

    densityEstimator = new CrowdingDistanceDensityEstimator<>();
    ranking = new FastNonDominatedSortRanking<>();
    replacement = new RankingAndDensityEstimatorReplacement<>(
        ranking, densityEstimator, Replacement.RemovalPolicy.ONE_SHOT);

    int tournamentSize = 2;
    selection = new NaryTournamentSelection<>(
        tournamentSize, variation.getMatingPoolSize(),
        new MultiComparator<>(
            Arrays.asList(
                Comparator.comparing(ranking::getRank),
                Comparator.comparing(densityEstimator::value).reversed())));

    createInitialPopulation = new PreComputedStructMSAsSolutionsCreation(problem, populationSize);

    if (numCores > 1) {
      evaluation = new MultiThreadedEvaluation<StructuralTM_MSASolution>(numCores, problem);
    } else {
      evaluation = new SequentialEvaluation<StructuralTM_MSASolution>(problem);
    }

    termination = new TerminationByEvaluations(maxEvaluations);
  }

  public StructuralTM_M2AlignBuilder setEvaluation(Evaluation<StructuralTM_MSASolution> evaluation) {
    this.evaluation = evaluation;

    return this;
  }

  public StructuralTM_M2Align build() {
    return new StructuralTM_M2Align(createInitialPopulation,
        evaluation, termination,
        selection, variation, replacement);
  }


 

}
