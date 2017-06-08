package studies.membrane.pdbtm;

import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionType;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import studies.membrane.MembraneConstants;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Statistics for transmembrane regions, grouped by interaction type.
 * Created by bittrich on 6/8/17.
 */
public class T04_StatisticsTransmembraneByInteractionType {
    public static void main(String[] args) {
        Map<Class<? extends PLIPInteraction>, StatisticsCollector.AminoAcidSummary> statistics = new HashMap<>();
        Stream.of(PLIPInteractionType.values())
                .forEach(plipInteractionType -> statistics.put(plipInteractionType.getDescribingClass(), new StatisticsCollector.AminoAcidSummary()));
        MembraneConstants.PdbtmAlphaNr.getInteractionsTransmembrane()
                .sequential()
                .forEach(plipInteraction -> {
                    StatisticsCollector.AminoAcidSummary aminoAcidSummary = statistics.get(plipInteraction.getClass());
                    aminoAcidSummary.accept((AminoAcid) plipInteraction.getPartner1());
                    aminoAcidSummary.accept((AminoAcid) plipInteraction.getPartner2());
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey().getSimpleName() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "type\t" + statistics.values().iterator().next().getHeaderLine() + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("transmembrane_by_interactionType.tsv"), output);
    }
}
