package design.aggregation.topology.gmlvq;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import design.ArffWriter;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static design.aggregation.S06_ExtractSequences.TOPOLOGIES;

/**
 * Steffen's approach to bin motifs merely by length and try to predict the topology based on that.
 * Created by S on 02.11.2016.
 */
@SuppressWarnings("unused")
public class S01_ConvertSequencesToArffByLength {
    // naive mapping of amino acids to their ordinal numbering of the amino acid enum
    private static final ToIntFunction<char[]> NAIVE_MAPPING = name -> determineAminoAcid(name).ordinal();
    // mapping to gutteridge grouping
    private static final ToIntFunction<char[]> GUTTERIDGE_MAPPING = name ->
            determineAminoAcid(name).getGutteridgeGrouping().ordinal();
    // mapping to membrane-preferring grouping
    private static final ToIntFunction<char[]> MEMBRANE_MAPPING = name ->
            determineAminoAcid(name).getANVILGrouping().ordinal();

    // switch between the way how amino acids are represented
    private static final ToIntFunction<char[]> MAPPING_RULE = NAIVE_MAPPING;

    public static void main(String[] args) throws IOException {
        Arrays.stream(SequenceMotifDefinition.values())
              // map to motif length
              .map(SequenceMotifDefinition::name)
              .map(motifName -> motifName.substring(2))
              .forEach(S01_ConvertSequencesToArffByLength::handleMotifLength);
    }

    private static AminoAcidFamily determineAminoAcid(char[] name) {
        return AminoAcidFamily.valueOfIgnoreCase(String.valueOf(name)).get();
    }

    private static void handleMotifLength(String lengthString) {
        // collect data
        List<String> data = TOPOLOGIES.stream()
            // skip transition regions for now
            .filter(topology -> topology.contains("tm"))
            .flatMap(topology -> handleTopology(lengthString, topology))
            .collect(Collectors.toList());

        ArffWriter.writeToArff(lengthString, data, Paths.get(DesignConstants.MOTIF_SEQUENCES_BY_LENGTH_DIR +
                lengthString + "-naive" + DesignConstants.ARFF_SUFFIX));
    }

    private static Stream<String> handleTopology(String lengthString, String topology) {
//        System.out.println("handling " + lengthString + " with " + topology);
        try {
            return Files.list(Paths.get(DesignConstants.EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR + topology + "/"))
                    // ensure correct fragment size
                    .filter(path -> path.getFileName().toString().substring(2).split("\\.")[0].equals(lengthString))
                    // read file
                    .flatMap(S01_ConvertSequencesToArffByLength::lines)
                    // convert to arff line
                    .map(line -> line.chars()
                                     .mapToObj(Character::toChars)
                                     .mapToInt(MAPPING_RULE)
                                     .mapToObj(String::valueOf)
                                     .collect(Collectors.joining(",", "", "," + topology)));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Stream<String> lines(Path path) {
        try {
            return Files.lines(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
