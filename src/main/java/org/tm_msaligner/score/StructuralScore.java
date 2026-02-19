package org.tm_msaligner.score;

import org.tm_msaligner.solution.StructuralTM_MSASolution;
import org.uma.jmetal.util.naming.DescribedEntity;
import org.tm_msaligner.util.AA;

public interface StructuralScore extends DescribedEntity {
  <S extends StructuralTM_MSASolution> double compute(S solution, AA[][]decodedSequences);
  boolean isAMinimizationScore();
  String getName();
}