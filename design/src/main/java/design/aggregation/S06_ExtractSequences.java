package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Extracts all sequence of 1 sequence motif for a given topology and writes them into 1 file.
 * Created by S on 02.11.2016.
 */
public class S06_ExtractSequences {
    public static final List<String> TOPOLOGIES = Arrays.asList("tm", "ntm", "trans");

    public static void main(String[] args) throws IOException {
        TOPOLOGIES.forEach(S06_ExtractSequences::handleTopology);
    }

    private static void handleTopology(String topology) {
        Arrays.stream(SequenceMotifDefinition.values())
            .forEach(definition -> {
                System.out.println("composing and writing " + topology + " with " + definition);
                String sequences = handleMotifDefinition(definition, topology);
                try {
                    Files.write(Paths.get(DesignConstants.EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR + topology + "/" +
                            definition + DesignConstants.SEQUENCE_SUFFIX), sequences.getBytes());
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                }
            });
    }

    private static String handleMotifDefinition(SequenceMotifDefinition definition, String topology) {
        try {
            return Files.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR))
                // skip dirs of false topology
                .filter(path -> path.getFileName().toString().startsWith(topology))
                // map to all files in dir
                .flatMap(S06_ExtractSequences::listFiles)
                // skip files of false motif defintion
                .filter(path -> path.getFileName().toString().startsWith(definition.name()))
                // parse fragment
                .map(ProteinParser::parsePDBFile)
                // extract sequence
                .map(Protein::getAminoAcidSequence)
                .collect(Collectors.joining(System.lineSeparator()));
        } catch (IOException e) {
           throw new UncheckedIOException(e);
        }
    }

    public static Stream<Path> listFiles(Path path) {
        try {
            return Files.list(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
