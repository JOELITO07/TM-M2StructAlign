package org.tm_msaligner.crossover;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.uma.jmetal.operator.crossover.CrossoverOperator;
import org.uma.jmetal.util.errorchecking.Check;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.tm_msaligner.solution.StructuralTM_MSASolution;


public class BioSPXMSACrossover implements CrossoverOperator<StructuralTM_MSASolution> {
  private final JMetalRandom randomGenerator;
  private final double probability;
  private final double alpha;

  public BioSPXMSACrossover(double probability , double alpha) {
    Check.probabilityIsValid(probability);
    Check.probabilityIsValid(alpha);
    
    this.randomGenerator = JMetalRandom.getInstance();
    this.probability = probability;
    this.alpha = alpha;

  }

  /**
   * Checks conditions and return the result of performing single point crossover
   *
   * @param parents
   * @return
   */
  public List<StructuralTM_MSASolution> execute(List<StructuralTM_MSASolution> parents) {
    Check.notNull(parents);
    Check.that(parents.size() == 2, "The number of parents is not 2");
    Check.that(
        parents.get(0).variables().size() == parents.get(1).variables().size(),
        "The two parents have different length: "
            + parents.get(0).variables().size()
            + ", "
            + parents.get(1).variables().size());

    return doCrossover(parents);
  }

  /**
   * Performs a single point crossover of two parents. Uses the same cutting point for all sequences
   *
   * @return a list containing the generated offspring
   */
  private List<StructuralTM_MSASolution> doCrossover(List<StructuralTM_MSASolution> parents) {
    StructuralTM_MSASolution parent1 = parents.get(0);
    StructuralTM_MSASolution parent2 = parents.get(1);

    List<StructuralTM_MSASolution> children = new ArrayList<StructuralTM_MSASolution>();

    children.add(MSACrossover(parent1, parent2));
    children.add(MSACrossover(parent2, parent1));

    return children;
  }

  private StructuralTM_MSASolution MSACrossover(StructuralTM_MSASolution parentA, StructuralTM_MSASolution parentB) {
    StructuralTM_MSASolution child;
    if (this.randomGenerator.nextDouble() < this.probability) {
      //int cut = selectRandomColumn(parentA);
      int cut = selectCutWithScanning(parentA);

      List<List<Integer>> gapsGroupFirstBloq = new ArrayList<List<Integer>>();
      List<Integer> carsCounterParentA = new ArrayList<Integer>();
      List<Integer> gapsGroup;
      int numgaps;

      for (int i = 0; i < parentA.variables().size(); i++) {
        gapsGroup = parentA.variables().get(i);

        List<Integer> gaps = new ArrayList<Integer>();
        numgaps = 0;
        for (int j = 1; j < gapsGroup.size(); j += 2) {
          if (cut >= gapsGroup.get(j)) {
            gaps.add(gapsGroup.get(j - 1));
            gaps.add(gapsGroup.get(j));
            numgaps += gapsGroup.get(j) - gapsGroup.get(j - 1) + 1;
          } else {
            if (cut >= gapsGroup.get(j - 1)) {
              gaps.add(gapsGroup.get(j - 1));
              gaps.add(cut);
              numgaps += cut - gapsGroup.get(j - 1) + 1;
            }
            break;
          }
        }
        gapsGroupFirstBloq.add(gaps);
        carsCounterParentA.add(cut - numgaps + 1);
      }

      int carsCountParentB;
      List<Integer> positionsToCutParentB = new ArrayList<Integer>();

      for (int i = 0; i < parentB.variables().size(); i++) {
        gapsGroup = parentB.variables().get(i);

        if (gapsGroup.size() > 0) {
          carsCountParentB = 0;
          for (int j = 0; j < gapsGroup.size(); j += 2) {
            if (j > 0) carsCountParentB += gapsGroup.get(j) - gapsGroup.get(j - 1) - 1;
            else carsCountParentB += gapsGroup.get(j);

            if (carsCountParentB >= carsCounterParentA.get(i)) {
              positionsToCutParentB.add(
                  gapsGroup.get(j) - (carsCountParentB - carsCounterParentA.get(i)));
              break;
            }
          }

          if (carsCountParentB < carsCounterParentA.get(i)) {
            if (gapsGroup.size() > 0) {
              carsCountParentB =
                  gapsGroup.get(gapsGroup.size() - 1)
                      + (carsCounterParentA.get(i) - carsCountParentB)
                      + 1;
              // if(carsCountParentB >= parent2.sizeAligment )
              // carsCountParentB=parent2.sizeAligment-1;
              positionsToCutParentB.add(carsCountParentB);
            }
          }
        } else { // SeqB has not Gaps
          positionsToCutParentB.add(carsCounterParentA.get(i));
        }
      }

      Integer MinPos = Collections.min(positionsToCutParentB);
      int pos;
      List<Integer> gaps;
      int lastGap, posA;
      for (int i = 0; i < parentB.variables().size(); i++) {
        posA = cut;
        pos = positionsToCutParentB.get(i);
        gaps = gapsGroupFirstBloq.get(i);
        if (pos > MinPos) {
          if (gaps.size() > 0) {
            lastGap = gaps.get(gaps.size() - 1);
            if (lastGap != posA) {
              gaps.add(posA + 1);
              gaps.add(posA + (pos - MinPos));
            } else {
              gaps.set(gaps.size() - 1, posA + (pos - MinPos));
            }
          } else {
            gaps.add(posA + 1);
            gaps.add(posA + (pos - MinPos));
          }

          posA += (pos - MinPos);
        }

        gapsGroup = parentB.variables().get(i);
        for (int j = 0; j < gapsGroup.size(); j += 2) {

          if (gapsGroup.get(j) >= pos) {
            gaps.add(posA + (gapsGroup.get(j) - pos) + 1);
            gaps.add(posA + (gapsGroup.get(j + 1) - pos) + 1);
          }
        }
      }

      child = new StructuralTM_MSASolution(parentA.getStructuralTMMSAProblem(), gapsGroupFirstBloq);

      child.mergeGapsGroups();

     

    } else {

      child = new StructuralTM_MSASolution(parentA);
    }

    
    return child;
  }

  /** Select a column randomly */
 /* public int selectRandomColumn(StructuralTM_MSASolution solution) {
    return randomGenerator.nextInt(1, solution.getAlignmentLength() - 1);
  }*/

  private int selectCutWithScanning(StructuralTM_MSASolution sol) {
        int L = sol.getAlignmentLength();
        int maxTries = 50;

        for (int t = 0; t < maxTries; t++) {
          int c = randomGenerator.nextInt(1, L - 1);

          if (sol.isGapColumn(c)) continue;
          boolean R = sol.riskIndicatorR(c);
          if (!R) return c;

          if (randomGenerator.nextDouble() < (1.0 - alpha)) return c;
        }

        return randomGenerator.nextInt(1, L - 1);
}

  public int numberOfRequiredParents() {
    return 2;
  }

  public int numberOfGeneratedChildren() {
    return 2;
  }

  public double crossoverProbability() {
	 return probability;
}


}
