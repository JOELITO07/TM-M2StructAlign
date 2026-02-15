package org.tm_msaligner.solution;


import java.util.List;
import org.tm_msaligner.problem.impl.MultiObjStructuralTMMSAProblem;
import org.tm_msaligner.util.AAArray;

public class StructuralTM_MSASolution extends TM_MSASolution {

    private final MultiObjStructuralTMMSAProblem problem ;

  public StructuralTM_MSASolution(MultiObjStructuralTMMSAProblem problem) {
    super(problem) ;
    this.problem = problem ;  
  }

  public StructuralTM_MSASolution(List<AAArray> AlignedSeqs, MultiObjStructuralTMMSAProblem problem) {
    super(AlignedSeqs, problem) ;
    this.problem = problem ;
  }

  public StructuralTM_MSASolution(MultiObjStructuralTMMSAProblem problem, List<List<Integer>> gapsGroups) {
    super(problem, gapsGroups) ;

    this.problem = problem ;

  }

  /** Copy Constructor */
  public StructuralTM_MSASolution(StructuralTM_MSASolution solution) {
    super(solution);
    problem = solution.problem ;

  }

   public float[][] getDistanceMatrixByIndex(int index) {     return problem.getDistanceMatrixByIndex(index);  }


}
