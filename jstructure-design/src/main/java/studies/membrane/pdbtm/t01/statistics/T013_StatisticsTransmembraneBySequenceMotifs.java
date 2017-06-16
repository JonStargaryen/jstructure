package studies.membrane.pdbtm.t01.statistics;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import studies.membrane.MembraneConstants;
import studies.membrane.pdbtm.StatisticsCollector;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Statistics for transmembrane regions, grouped by sequence motifs.
 * Created by bittrich on 6/8/17.
 */
public class T013_StatisticsTransmembraneBySequenceMotifs {
    public static void main(String[] args) {
        Map<SequenceMotifDefinition, StatisticsCollector.AminoAcidSummary> statistics = new HashMap<>();
        Stream.of(SequenceMotifDefinition.values())
                .forEach(sequenceMotifDefinition -> statistics.put(sequenceMotifDefinition, new StatisticsCollector.AminoAcidSummary()));
        MembraneConstants.PdbtmAlphaNr.getSequenceMotifsTransmembrane()
                .sequential()
                .forEach(sequenceMotif -> {
                    StatisticsCollector.AminoAcidSummary aminoAcidSummary = statistics.get(sequenceMotif.getMotifDefinition());
                    sequenceMotif.getAminoAcids().forEach(aminoAcidSummary);
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey().name() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "motif\t" + statistics.values().iterator().next().getHeaderLine() + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("transmembrane_by_sequenceMotif.tsv"), output);
    }
}
