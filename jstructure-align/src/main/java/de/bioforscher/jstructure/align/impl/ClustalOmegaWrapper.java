package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.align.MultipleSequenceAlignmentResult;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Employs a local clustal omega executable to run the alignment.
 * Created by bittrich on 7/14/17.
 */
public class ClustalOmegaWrapper extends AbstractClustalOmega {
    private static final String CLUSTALOMEGA_COMMAND = "clustalo";

    @Override
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
}
