package studies.membrane.pdbtm.t01;

import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionType;
import studies.membrane.MembraneConstants;
import studies.membrane.pdbtm.StatisticsCollector;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Extract interactions between amino acids in different helices, grouped by interaction type.
 * Created by bittrich on 6/8/17.
 */
public class T017_StatisticsForInteractingHelicesByInteractionType {
    public static void main(String[] args) {
        Map<Class<? extends PLIPInteraction>, StatisticsCollector.InteractingAminoAcidSummary> statistics = new HashMap<>();
        Stream.of(PLIPInteractionType.values())
                .forEach(plipInteractionType -> statistics.put(plipInteractionType.getDescribingClass(), new StatisticsCollector.InteractingAminoAcidSummary()));
        MembraneConstants.PdbtmAlphaNr.getInteractionsBetweenTransmembraneHelices()
                .sequential()
                .forEach(plipInteraction -> {
                    StatisticsCollector.InteractingAminoAcidSummary aminoAcidSummary = statistics.get(plipInteraction.getClass());
                    aminoAcidSummary.accept(plipInteraction);
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey().getSimpleName() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "type\t" + statistics.values().iterator().next().getHeaderLine() + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("interactingHelices_by_interactionType.tsv"), output);
    }
}
