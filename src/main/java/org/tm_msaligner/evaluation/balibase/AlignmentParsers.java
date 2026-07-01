package org.tm_msaligner.evaluation.balibase;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Parsers for BAliBASE XML and common MSA text formats. */
public final class AlignmentParsers {
    private AlignmentParsers() {
    }

    public static Alignment readReference(Path path) throws IOException {
        String content = Files.readString(path, StandardCharsets.UTF_8);
        return looksLikeXml(content) ? readXml(content) : readAlignment(content);
    }

    public static Alignment readTest(Path path) throws IOException {
        return readAlignment(Files.readString(path, StandardCharsets.UTF_8));
    }

    public static Alignment readAlignment(String content) throws IOException {
        String trimmed = content.stripLeading();
        if (looksLikeXml(trimmed)) {
            return readXml(trimmed);
        }
        if (trimmed.startsWith(">")) {
            return readFasta(trimmed);
        }
        if (trimmed.toUpperCase(Locale.ROOT).startsWith("CLUSTAL")) {
            return readClustal(trimmed);
        }
        return readMsf(trimmed);
    }

    private static boolean looksLikeXml(String content) {
        String text = content.stripLeading();
        return text.startsWith("<?xml") || text.startsWith("<ms_alignment") || text.startsWith("<alignment");
    }

    /** Reads BAliBASE XML and extracts seq-name, seq-data, and optional coreblock column-score. */
    public static Alignment readXml(String content) throws IOException {
        try {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            factory.setFeature(XMLConstants.FEATURE_SECURE_PROCESSING, true);
            disableFeature(factory, "http://apache.org/xml/features/disallow-doctype-decl");
            factory.setExpandEntityReferences(false);
            Document document = factory.newDocumentBuilder().parse(new InputSource(new StringReader(content)));
            document.getDocumentElement().normalize();

            List<AlignmentSequence> sequences = new ArrayList<>();
            NodeList sequenceNodes = document.getElementsByTagName("sequence");
            for (int i = 0; i < sequenceNodes.getLength(); i++) {
                Element sequence = (Element) sequenceNodes.item(i);
                String name = firstText(sequence, "seq-name");
                String data = firstText(sequence, "seq-data");
                if (name != null && data != null) {
                    sequences.add(new AlignmentSequence(name, data));
                }
            }
            if (sequences.isEmpty()) {
                throw new IOException("No <sequence> entries with <seq-name> and <seq-data> were found in XML reference");
            }
            return new Alignment(sequences, readCoreBlockMask(document));
        } catch (Exception e) {
            if (e instanceof IOException) {
                throw (IOException) e;
            }
            throw new IOException("Cannot parse BAliBASE XML reference", e);
        }
    }

    private static void disableFeature(DocumentBuilderFactory factory, String feature) {
        try {
            factory.setFeature(feature, true);
        } catch (Exception ignored) {
            // Some JDK XML providers do not expose every hardening feature.
        }
    }

    private static String firstText(Element parent, String tagName) {
        NodeList nodes = parent.getElementsByTagName(tagName);
        if (nodes.getLength() == 0) {
            return null;
        }
        String text = nodes.item(0).getTextContent();
        return text == null ? null : text.trim();
    }

    private static int[] readCoreBlockMask(Document document) {
        NodeList scoreNodes = document.getElementsByTagName("column-score");
        for (int i = 0; i < scoreNodes.getLength(); i++) {
            Element score = (Element) scoreNodes.item(i);
            String name = firstText(score, "colsco-name");
            String data = firstText(score, "colsco-data");
            if (name != null && data != null && "coreblock".equalsIgnoreCase(name.trim())) {
                return Arrays.stream(data.trim().split("\\s+"))
                        .filter(token -> !token.isBlank())
                        .mapToInt(Integer::parseInt)
                        .toArray();
            }
        }
        return null;
    }

    /** Reads GCG/MSF format, the native format expected by the original bali_score test alignment. */
    public static Alignment readMsf(String content) throws IOException {
        int marker = content.indexOf("//");
        if (marker < 0) {
            throw new IOException("Input does not look like MSF/GCG format: missing // alignment marker");
        }
        Map<String, StringBuilder> sequences = new LinkedHashMap<>();
        try (BufferedReader reader = new BufferedReader(new StringReader(content.substring(marker + 2)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                String trimmed = line.trim();
                if (trimmed.isEmpty() || Character.isDigit(trimmed.charAt(0))) {
                    continue;
                }
                String[] parts = trimmed.split("\\s+", 2);
                if (parts.length < 2 || !containsResidues(parts[1])) {
                    continue;
                }
                String name = parts[0];
                String residues = parts[1].replaceAll("\\s+", "");
                sequences.computeIfAbsent(name, ignored -> new StringBuilder()).append(residues);
            }
        }
        return toAlignment(sequences, "MSF/GCG");
    }

    public static Alignment readFasta(String content) throws IOException {
        Map<String, StringBuilder> sequences = new LinkedHashMap<>();
        String currentName = null;
        try (Reader sr = new StringReader(content); BufferedReader reader = new BufferedReader(sr)) {
            String line;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }
                if (trimmed.startsWith(">")) {
                    currentName = trimmed.substring(1).trim().split("\\s+", 2)[0];
                    sequences.putIfAbsent(currentName, new StringBuilder());
                } else if (currentName != null) {
                    sequences.get(currentName).append(trimmed);
                }
            }
        }
        return toAlignment(sequences, "FASTA");
    }

    public static Alignment readClustal(String content) throws IOException {
        Map<String, StringBuilder> sequences = new LinkedHashMap<>();
        try (BufferedReader reader = new BufferedReader(new StringReader(content))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank() || line.toUpperCase(Locale.ROOT).startsWith("CLUSTAL")) {
                    continue;
                }
                if (Character.isWhitespace(line.charAt(0)) || line.trim().matches("^[*. :]+$")) {
                    continue;
                }
                String[] parts = line.trim().split("\\s+", 3);
                if (parts.length >= 2 && containsResidues(parts[1])) {
                    sequences.computeIfAbsent(parts[0], ignored -> new StringBuilder()).append(parts[1]);
                }
            }
        }
        return toAlignment(sequences, "CLUSTAL");
    }

    private static boolean containsResidues(String text) {
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (Character.isLetter(c) || c == '-' || c == '.' || c == '~' || c == '*') {
                return true;
            }
        }
        return false;
    }

    private static Alignment toAlignment(Map<String, StringBuilder> sequences, String formatName) throws IOException {
        if (sequences.isEmpty()) {
            throw new IOException("No sequences were parsed from " + formatName + " alignment");
        }
        List<AlignmentSequence> list = new ArrayList<>(sequences.size());
        for (Map.Entry<String, StringBuilder> entry : sequences.entrySet()) {
            list.add(new AlignmentSequence(entry.getKey(), entry.getValue().toString()));
        }
        return new Alignment(list, null);
    }
}
