package org.tm_msaligner.solutionscreation;

import java.sql.Struct;
import java.util.ArrayList;
import java.util.List;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.tm_msaligner.crossover.BioSPXMSACrossover;
import org.tm_msaligner.crossover.SPXMSACrossover;
import org.tm_msaligner.problem.StandardTMMSAProblem;
import org.tm_msaligner.problem.StructuralTMMSAProblem;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.solution.TM_MSASolution;
import org.tm_msaligner.util.AAArray;

public class PreComputedStructMSAsSolutionsCreation implements SolutionsCreation<StructuralTM_MSASolution> {

  private final int numberOfSolutionsToCreate;
  private final StructuralTMMSAProblem problem;

  public PreComputedStructMSAsSolutionsCreation(StructuralTMMSAProblem problem,
      int numberOfSolutionsToCreate) {
    this.numberOfSolutionsToCreate = numberOfSolutionsToCreate;
    this.problem = problem;
  }

  @Override
  public List<StructuralTM_MSASolution> create() {
    List<StructuralTM_MSASolution> population = new ArrayList<>(numberOfSolutionsToCreate);
    JMetalRandom randomGenerator = JMetalRandom.getInstance();

    for (List<AAArray> sequenceList : problem.listOfPrecomputedStringAlignments) {
      StructuralTM_MSASolution newIndividual = new StructuralTM_MSASolution(sequenceList, problem);
      population.add(newIndividual);
    }

    int parent1, parent2;
    List<StructuralTM_MSASolution> children, parents;
    BioSPXMSACrossover crossover = new BioSPXMSACrossover(1);

    while (population.size() < numberOfSolutionsToCreate) {
      parents = new ArrayList<StructuralTM_MSASolution>();

      parent1 = randomGenerator.nextInt(0, population.size() - 1);
      do {
        parent2 = randomGenerator.nextInt(0, population.size() - 1);
      } while (parent1 == parent2);
      parents.add(population.get(parent1));
      parents.add(population.get(parent2));

      children = crossover.execute(parents);

      population.add(children.get(0));
      population.add(children.get(1));
    }

    return population;
  }
}
