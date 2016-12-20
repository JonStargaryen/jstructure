package design.learning;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static design.DesignConstants.DELIMITER;

/**
 * Converts a cluster summary  file to its *.arff representation.
 * Created by bittrich on 12/19/16.
 */
public class ClusterSummaryToArffConverter {
    private static final String NOMINAL_AMINO_ACID_STRING = Arrays.stream(AminoAcid.values())
            .map(AminoAcid::getOneLetterCode)
            .collect(Collectors.joining(",", "{", "}"));

    // naive mapping of amino acids to their ordinal numbering of the amino acid enum
    private static final ToIntFunction<char[]> NAIVE_MAPPING = name -> determineAminoAcid(name).ordinal();

    private static AminoAcid determineAminoAcid(char[] name) {
        return AminoAcid.valueOfIgnoreCase(String.valueOf(name));
    }

    /**
     * Exports the given data to a String. The last entry in the data argument will serve as class label.
     * @param inputDataPath the data container
     * @return the String representation of the data
     */
    public static String convertToArff(Path inputDataPath) {
        try {
            // read file
            List<String> lines = Files.readAllLines(inputDataPath);
            // drop first line
            lines.remove(0);

            // extract classes
            Set<String> uniqueClassLabels = lines.stream()
                    .map(line -> line.split(DELIMITER))
                    .map(split -> split[0])
                    .collect(Collectors.toSet());

            String sequence = lines.stream()
                    .limit(1)
                    .map(line -> line.split(DELIMITER))
                    .map(split -> split[4])
                    .findFirst()
                    .orElseThrow(NoSuchElementException::new);

            int numberOfVariablePositions = sequence.length() - 2;
            String motifString = sequence.charAt(0) + "" + sequence.charAt(sequence.length() - 1) + (numberOfVariablePositions + 1);

            return lines.stream()
                    .map(line -> line.split(DELIMITER))
                    .map(split -> convertToNumericSequenceRepresentation(split[4]) + "," + split[0])
                    .collect(Collectors.joining(System.lineSeparator(),
                            "@RELATION " + motifString + System.lineSeparator() +
                                    IntStream.range(0, numberOfVariablePositions)
                                            .mapToObj(index -> "@ATTRIBUTE position" + index + " NUMERIC")
                                            .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                                    "@ATTRIBUTE class " + uniqueClassLabels.stream()
                                            .collect(Collectors.joining(",",
                                                    "{",
                                                    "}")) + System.lineSeparator() +
                                    "@DATA" + System.lineSeparator(),
                            ""));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String convertToNumericSequenceRepresentation(String sequence) {
        // trim to variable positions
        sequence = sequence.substring(1, sequence.length() - 1);

        return sequence.chars()
                .mapToObj(Character::toChars)
                .mapToInt(NAIVE_MAPPING)
                .mapToObj(String::valueOf)
                .collect(Collectors.joining(","));
    }

    /**
     * Exports the given data to the file system.
     * @param inputDataPath the data container - each entry must provide a suitable {@link Object#toString()} implementation
     * @param outputDataPath where to write?
     * @see #convertToArff(Path)
     */
    public static void writeToArff(Path inputDataPath, Path outputDataPath) {
        DesignConstants.write(outputDataPath, convertToArff(inputDataPath).getBytes());
    }
}
