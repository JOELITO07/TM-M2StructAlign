package org.tm_msaligner.evaluation;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Locale;
import org.tm_msaligner.util.AAArray;

/** Reads the project's name/sequence/topology three-line format. */
final class TopologyDataset {
  private final LinkedHashMap<String, AAArray> residuesByProtein;

  private TopologyDataset(LinkedHashMap<String, AAArray> residuesByProtein) {
    this.residuesByProtein = residuesByProtein;
  }

  static TopologyDataset read(Path path) throws IOException {
    LinkedHashMap<String, AAArray> proteins = new LinkedHashMap<>();

    try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
      while (true) {
        String header = nextNonEmptyLine(reader);
        if (header == null) {
          break;
        }
        if (!header.startsWith(">")) {
          throw new IOException(
              "Expected a FASTA header in topology file " + path + ", obtained: " + header);
        }

        String sequence = nextNonEmptyLine(reader);
        String topology = nextNonEmptyLine(reader);
        if (sequence == null || topology == null) {
          throw new IOException("Incomplete three-line record for " + header + " in " + path);
        }

        sequence = removeWhitespace(sequence).toUpperCase(Locale.ROOT);
        topology = removeWhitespace(topology).toUpperCase(Locale.ROOT);
        if (sequence.length() != topology.length()) {
          throw new IOException(
              "Sequence/topology length mismatch for " + header + " in " + path
                  + ": sequence=" + sequence.length() + ", topology=" + topology.length());
        }
        if (sequence.indexOf('-') >= 0 || sequence.indexOf('.') >= 0) {
          throw new IOException(
              "Topology dataset sequences must be ungapped; found a gap in " + header);
        }

        String id = FastaAlignment.normalizeProteinId(header);
        AAArray residues = new AAArray(sequence, topology);
        for (int index = 0; index < sequence.length(); index++) {
          residues.AAAt(index).setSeqIndex(index);
        }
        if (proteins.putIfAbsent(id, residues) != null) {
          throw new IOException(
              "Duplicated normalized protein ID '" + id + "' in topology file " + path);
        }
      }
    }

    if (proteins.isEmpty()) {
      throw new IOException("No topology records found in " + path);
    }
    return new TopologyDataset(proteins);
  }

  private static String nextNonEmptyLine(BufferedReader reader) throws IOException {
    String line;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      if (!line.isEmpty()) {
        return line;
      }
    }
    return null;
  }

  private static String removeWhitespace(String value) {
    return value.replaceAll("\\s+", "");
  }

  AAArray residues(String proteinId) {
    return residuesByProtein.get(proteinId);
  }
}
