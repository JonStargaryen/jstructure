package studies.membrane.pdbtm;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import studies.membrane.MembraneConstants;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Statistics for transmembrane regions, grouped by sequence motifs.
 * Created by bittrich on 6/8/17.
 */
public class T03_StatisticsTransmembraneBySequenceMotifs {
    public static void main(String[] args) {
        Map<SequenceMotifDefinition, StatisticsCollector.OccurrenceSummary> statistics = new HashMap<>();
        Stream.of(SequenceMotifDefinition.values())
                .forEach(sequenceMotifDefinition -> statistics.put(sequenceMotifDefinition, new StatisticsCollector.OccurrenceSummary()));
        MembraneConstants.getSequenceMotifsOfPdbtmAlphaNrListTransmembrane()
                .sequential()
                .forEach(sequenceMotif -> {
                    StatisticsCollector.OccurrenceSummary occurrenceSummary = statistics.get(sequenceMotif.getMotifDefinition());
                    sequenceMotif.getAminoAcids().forEach(occurrenceSummary);
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey().name() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "motif\t" + StatisticsCollector.OccurrenceSummary.getHeaderLine(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("transmembrane_by_sequenceMotif.tsv"), output);
    }
}
