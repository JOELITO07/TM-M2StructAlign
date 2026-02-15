package org.tm_msaligner.problem;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.tm_msaligner.solution.TM_MSASolution;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.errorchecking.JMetalException;

import org.tm_msaligner.util.AAArray;

public class StandardTMMSAProblem extends AbstractGenericTM_MSAProblem<TM_MSASolution> {

  public List<AAArray> originalSequences;
  public List<List<AAArray>> listOfPrecomputedStringAlignments;
  public List<StringBuilder> listOfSequenceNames;
  public long[] MaxMinSegmentAlignScore;

  public StandardTMMSAProblem(String msaProblemFileName, List<String> preComputedFiles)
      throws IOException {

    if (preComputedFiles.size() < 2) {
      throw new JMetalException(
          "Wrong number of Pre-computed Alignments, Minimum 2 files are required");
    }

    readSequenceFromFile(msaProblemFileName);
    setNumberOfVariables(originalSequences.size());
    MaxMinSegmentAlignScore = getMaxMinScoreSegmentAlign();
    listOfPrecomputedStringAlignments = readPreComputedAlignments(preComputedFiles);
  }

   public List<List<AAArray>> readPreComputedAlignments(List<String> dataFiles) {
    List<List<AAArray>> listPreAlignments = new ArrayList<List<AAArray>>();
    for (String dataFile : dataFiles) {
      try {
        List<AAArray> seqAligned = readDataFromFastaFile(dataFile);
        TM_MSASolution sol = new TM_MSASolution(seqAligned, this);
        if (!sol.isValid()) {
          JMetalLogger.logger.warning("MSA in file " + dataFile + " is not Valid");
        } else {
          listPreAlignments.add(seqAligned);
        }
      } catch (Exception e) {
        throw new JMetalException(
            "Error reading data from fasta files " + dataFile + ". Message: " + e);
      }
    }

    if (listPreAlignments.size() < 2) {
      throw new JMetalException("More than one PreComputedAlignment is needed");
    }

    return listPreAlignments;
  }

  public List<AAArray> readDataFromFastaFile(String dataFile)
      throws IOException {

    List<AAArray> sequenceList = new ArrayList<AAArray>();

    LinkedHashMap<String, ProteinSequence> sequences =
        FastaReaderHelper.readFastaProteinSequence(new File(dataFile));

    for (Map.Entry<String, ProteinSequence> entry : sequences.entrySet()) {
      sequenceList.add(new AAArray(entry.getValue().getSequenceAsString()));
    }

    return sequenceList;
  }

  public List<StringBuilder> getListOfSequenceNames() {
    return listOfSequenceNames;
  }

  public List<StringBuilder> readSeqNameFromAlignment(String dataFile)
      throws IOException, CompoundNotFoundException {

    List<StringBuilder> SeqNameList = new ArrayList<StringBuilder>();

    LinkedHashMap<String, ProteinSequence> sequences =
        FastaReaderHelper.readFastaProteinSequence(new File(dataFile));

    for (Map.Entry<String, ProteinSequence> entry : sequences.entrySet()) {
      SeqNameList.add(new StringBuilder(entry.getValue().getOriginalHeader()));
    }

    return SeqNameList;
  }

  public void printSequenceListToFasta(List<AAArray> solutionList, String fileName)
      throws Exception {
    List<ProteinSequence> proteinSequences = new ArrayList<ProteinSequence>();
    for (AAArray sequence : solutionList) {
      proteinSequences.add(new ProteinSequence(sequence.toString()));
    }

    FastaWriterHelper.writeProteinSequence(new File(fileName), proteinSequences);
  }

  public long[] getMaxMinScoreSegmentAlign() {
    long[] MaxMinScores = new long[2];
    long MaxScore = 0, MinScore = 0;
    AAArray seq;
    int numSeqs = originalSequences.size();

    for (int i = 0; i < numSeqs - 1; i++) {
      seq = originalSequences.get(i);
      for (int l = 0; l < seq.getSize(); l++) {
        MinScore += (-1) * (numSeqs - i - 1);
        if (seq.AAAt(l).getType().isTMRegion()) {
          MaxScore += 4 * (numSeqs - i - 1);
        } else {
          MaxScore += 2 * (numSeqs - i - 1);
        }
      }
    }
    MaxMinScores[0] = MaxScore;
    MaxMinScores[1] = MinScore;
    return MaxMinScores;
  }

  public TM_MSASolution evaluate(TM_MSASolution tm_msaSolution) { return null; }
  public TM_MSASolution createSolution() { return null; }

  void readSequenceFromFile(String file) {
    originalSequences = new ArrayList<AAArray>();
    listOfSequenceNames = new ArrayList<StringBuilder>();

    try {
      BufferedReader in = new BufferedReader(new FileReader(file));

      int status = 0;
      String line, regiones;
      int posSepName;

      for (line = in.readLine(); line != null; line = in.readLine()) {
        line = line.trim();
        if (status == 0) {
          if (line.length() > 0 && line.charAt(0) == '>') {
            posSepName = line.indexOf('|');
            listOfSequenceNames.add(new StringBuilder(
                posSepName > 0 ? line.substring(1, posSepName) : line.substring(1)));
          } else {
            throw new IOException("Name of Sequence must starts with '>'");
          }
          status = 1;
        } else if (status == 1) {
          regiones = in.readLine().trim();
          if (regiones == null) {
            throw new IOException("Regions of Sequence is empty");
          }
          if (regiones.length() != line.length()) {
            throw new IOException("Regions of Sequence is empty");
          }

          originalSequences.add(new AAArray(line, regiones));
          status = 0;
        }
      }

      if (originalSequences.size() != listOfSequenceNames.size()) {
        throw new IOException("Names wiht Sequences are not equals");
      }

    } catch (IOException e) {
      System.out.println("Error when reading " + file);
      e.printStackTrace();
    }
  }

  public int getSizeOfOriginalSequence(int i) {
    if (i < originalSequences.size()) {
      return originalSequences.get(i).getSize();
    } else {
      System.out.println("Error getting size of Original Sequence " + i +
          " and the Number of Sequences is " + originalSequences.size());
      return 0;
    }
  }
}
