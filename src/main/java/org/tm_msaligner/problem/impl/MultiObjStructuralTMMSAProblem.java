package org.tm_msaligner.problem.impl;


import java.io.IOException;
import java.util.List;
import org.tm_msaligner.problem.StructuralTMMSAProblem;
import org.tm_msaligner.score.StructuralScore;
import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.tm_msaligner.util.AA;

public class MultiObjStructuralTMMSAProblem extends StructuralTMMSAProblem {

  private final List<StructuralScore> scoreList;

  public MultiObjStructuralTMMSAProblem(String msaProblemFileName, List<StructuralScore> scoreList,
      List<String> preComputedFiles, String distanceDatasetDir, String Name) throws IOException {
    super(msaProblemFileName, preComputedFiles, distanceDatasetDir);

    setNumberOfObjectives(scoreList.size());
    setName(Name); // "Multi Objective Structural TM-MSA Problem"

    this.scoreList = scoreList;
  }

  @Override
  public StructuralTM_MSASolution evaluate(StructuralTM_MSASolution solution) {
    solution.removeGapColumns();
    AA[][] decodedSequences = solution.decodeToMatrix();

    for (int i = 0; i < numberOfObjectives(); i++) {
      solution.objectives()[i] = scoreList.get(i).compute(solution, decodedSequences) *
          (scoreList.get(i).isAMinimizationScore() ? 1.0 : -1.0);
    }

    return solution;
  }
}
