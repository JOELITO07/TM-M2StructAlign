package org.tm_msaligner.evaluation;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** A small strict FASTA alignment reader that preserves gap characters. */
final class FastaAlignment {
  private static final Set<String> STRUCTURE_EXTENSIONS =
      Set.of(".pdb", ".ent", ".cif", ".mmcif");

  private final LinkedHashMap<String, String> sequences;
  private final int alignmentLength;

  private FastaAlignment(LinkedHashMap<String, String> sequences, int alignmentLength) {
    this.sequences = sequences;
    this.alignmentLength = alignmentLength;
  }

  static FastaAlignment read(Path path) throws IOException {
    LinkedHashMap<String, StringBuilder> parsed = new LinkedHashMap<>();
    String currentId = null;

    try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
      String line;
      int lineNumber = 0;
      while ((line = reader.readLine()) != null) {
        lineNumber++;
        line = line.trim();
        if (line.isEmpty() || line.startsWith(";")) {
          continue;
        }

        if (line.startsWith(">")) {
          currentId = normalizeProteinId(line);
          if (parsed.putIfAbsent(currentId, new StringBuilder()) != null) {
            throw new IOException(
                "Duplicated normalized protein ID '" + currentId + "' in " + path);
          }
        } else {
          if (currentId == null) {
            throw new IOException(
                "Sequence data before the first FASTA header in " + path
                    + " at line " + lineNumber);
          }
          parsed.get(currentId).append(removeWhitespace(line).replace('.', '-').toUpperCase(
              Locale.ROOT));
        }
      }
    }

    if (parsed.isEmpty()) {
      throw new IOException("No FASTA records found in " + path);
    }

    LinkedHashMap<String, String> sequences = new LinkedHashMap<>();
    int expectedLength = -1;
    for (Map.Entry<String, StringBuilder> entry : parsed.entrySet()) {
      String sequence = entry.getValue().toString();
      if (sequence.isEmpty()) {
        throw new IOException(
            "Empty sequence for protein '" + entry.getKey() + "' in " + path);
      }
      if (expectedLength < 0) {
        expectedLength = sequence.length();
      } else if (sequence.length() != expectedLength) {
        throw new IOException(
            "The FASTA records in " + path + " do not have a common alignment length: "
                + entry.getKey() + " has " + sequence.length()
                + " columns; expected " + expectedLength);
      }
      sequences.put(entry.getKey(), sequence);
    }

    return new FastaAlignment(sequences, expectedLength);
  }

  static String normalizeProteinId(String header) throws IOException {
    if (header == null) {
      throw new IOException("Null FASTA header");
    }

    String value = header.trim();
    if (value.startsWith(">")) {
      value = value.substring(1).trim();
    }
    if (value.isEmpty()) {
      throw new IOException("Empty FASTA identifier");
    }

    value = value.split("\\s+", 2)[0];
    String[] fields = value.split("\\|");
    if (fields.length >= 2
        && ("sp".equalsIgnoreCase(fields[0]) || "tr".equalsIgnoreCase(fields[0]))) {
      value = fields[1];
    } else {
      value = fields[0];
    }

    int chainSeparator = value.indexOf(':');
    if (chainSeparator >= 0) {
      value = value.substring(0, chainSeparator);
    }

    value = value.trim();
    String lowerValue = value.toLowerCase(Locale.ROOT);
    for (String extension : STRUCTURE_EXTENSIONS) {
      if (lowerValue.endsWith(extension)) {
        value = value.substring(0, value.length() - extension.length());
        break;
      }
    }

    value = value.trim().toLowerCase(Locale.ROOT);
    if (value.isEmpty()) {
      throw new IOException("Empty normalized FASTA identifier derived from: " + header);
    }
    return value;
  }

  private static String removeWhitespace(String value) {
    return value.replaceAll("\\s+", "");
  }

  Map<String, String> sequences() {
    return Collections.unmodifiableMap(sequences);
  }

  int alignmentLength() {
    return alignmentLength;
  }
}
