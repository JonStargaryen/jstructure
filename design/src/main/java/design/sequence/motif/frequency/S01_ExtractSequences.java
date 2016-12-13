package design.sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import design.ProteinSource;

import java.util.Map;
import java.util.stream.Collectors;

/**
 * Returns all sequence clusters (described the actual structural clusters). Is there a correlation between sequences
 * and multiple structure clusters?
 * Created by S on 13.12.2016.
 */
public class S01_ExtractSequences {
    public static void main(String[] args) {
        extractSequences("trans", SequenceMotifDefinition.VL4).entrySet().forEach(System.out::println);
    }

    static Map<String, String> extractSequences(String topology, SequenceMotifDefinition sequenceMotifDefinition) {
        return ProteinSource.getGroupedFragments(topology, sequenceMotifDefinition)
                .entrySet()
                .stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().stream()
                        // cast is safe as each fragment contains multiple groups
                        .map(GroupContainer.class::cast)
                        .map(GroupContainer::getAminoAcidSequence)
                        .collect(Collectors.joining(System.lineSeparator()))));
    }
}
