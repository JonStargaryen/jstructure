package design.sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.util.Collection;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static design.DesignConstants.DECIMAL_FORMAT;

/**
 * Converts sequences of clusters to their JSON representation to visualize them by sequence boxes.
 * Created by S on 13.12.2016.
 */
public class S03_ExtractedSequencesToJSON {
    public static void main(String[] args) {
        System.out.println(sequencesToJSON("trans", SequenceMotifDefinition.VL4));
    }

    private static Map<String, String> sequencesToJSON(String topology, SequenceMotifDefinition sequenceMotifDefinition) {
        Map<String, String> sequenceClusters = S01_ExtractSequences.extractSequences(topology, sequenceMotifDefinition);
        sequenceClusters.put(DesignConstants.RARE_CLUSTER_NAME, S02_ExtractMissingSequences.extractRareSequences(topology, sequenceMotifDefinition));

        return sequenceClusters.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry ->
                        toJSON(new SequenceMotifRepresentation(Pattern.compile("\n").splitAsStream(entry.getValue())
                                .peek(System.out::println)
                                .collect(Collectors.toList())))));
    }

    private static String toJSON(SequenceMotifRepresentation representation) {
        return representation.getFrequencies().stream()
                .map(Map::entrySet)
                .flatMap(Collection::stream)
                .map(Map.Entry::getValue)
                .map(DECIMAL_FORMAT::format)
                .collect(Collectors.joining(",",
                        System.lineSeparator() + "[",
                        "]" + System.lineSeparator()));
    }
}
