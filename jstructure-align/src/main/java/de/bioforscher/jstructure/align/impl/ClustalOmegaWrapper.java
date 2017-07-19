package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Employs a local clustal omega executable to run the alignment.
 * Created by bittrich on 7/14/17.
 */
public class ClustalOmegaWrapper implements MultipleSequenceAligner {
    private static final Pattern SEQUENCE_PATTERN = Pattern.compile(">");
    private static final Pattern LINE_PATTERN = Pattern.compile("\\s+");
    private static final String CLUSTALOMEGA_COMMAND = "clustalo";

    public MultipleSequenceAlignmentResult align(List<String> fastaSequences) throws AlignmentException {
        try {
            Path tmpDirectory = Files.createTempDirectory("clustalomega");

            // write input sequence file
            Path sequencePath = tmpDirectory.resolve("sequences.fasta");
            Files.write(sequencePath, fastaSequences.stream()
                    .collect(Collectors.joining(System.lineSeparator())).getBytes());

            // expected output file
            Path outputPath = tmpDirectory.resolve("output.fasta");

            ProcessBuilder processBuilder = new ProcessBuilder(CLUSTALOMEGA_COMMAND,
                    "-i",
                    sequencePath.toFile().getAbsolutePath(),
                    "-o",
                    outputPath.toFile().getAbsolutePath());
            processBuilder.start().waitFor();

            return createMultipleSequenceAlignmentResult(Files.lines(outputPath)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            throw new AlignmentException(e);
        }
    }

    private MultipleSequenceAlignmentResult createMultipleSequenceAlignmentResult(String rawAlignmentString) {
        // use LinkedHashMap to sort entries according their addition to the map
        Map<String, String> alignment = new LinkedHashMap<>();
        SEQUENCE_PATTERN.splitAsStream(rawAlignmentString)
                // skip first (empty) pair
                .skip(1)
                .forEach(line -> {
                    // split alignment at newlines
                    String[] split = LINE_PATTERN.split(line);
                    // substring to drop FASTA-character
                    alignment.put(split[0],
                            // skip id and join all other lines
                            Stream.of(split).skip(1).collect(Collectors.joining()));
                });

        return new MultipleSequenceAlignmentResultImpl(alignment);
    }
}
