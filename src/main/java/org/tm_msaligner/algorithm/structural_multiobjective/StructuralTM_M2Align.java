package org.tm_msaligner.algorithm.structural_multiobjective;

import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;
import org.uma.jmetal.component.catalogue.common.evaluation.Evaluation;
import org.uma.jmetal.component.catalogue.common.solutionscreation.SolutionsCreation;
import org.uma.jmetal.component.catalogue.common.termination.Termination;
import org.uma.jmetal.component.catalogue.ea.replacement.Replacement;
import org.uma.jmetal.component.catalogue.ea.selection.Selection;
import org.uma.jmetal.component.catalogue.ea.variation.Variation;

public class StructuralTM_M2Align extends EvolutionaryAlgorithm<StructuralTM_MSASolution> {

  public StructuralTM_M2Align(SolutionsCreation<StructuralTM_MSASolution> initialPopulationCreation,
      Evaluation<StructuralTM_MSASolution> evaluation,
      Termination termination,
      Selection<StructuralTM_MSASolution> selection,
      Variation<StructuralTM_MSASolution> variation,
      Replacement<StructuralTM_MSASolution> replacement) {
    super("StructuralTM-M2Align", initialPopulationCreation, evaluation, termination, selection, variation,replacement);

  }

}
