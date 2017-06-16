package studies.membrane.pdbtm.t01.statistics;

import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import studies.membrane.MembraneConstants;
import studies.membrane.pdbtm.StatisticsCollector;

import java.util.Map;
import java.util.stream.Collectors;

/**
 * Statistics now with a split by topology.
 * Created by bittrich on 6/7/17.
 */
public class T012_StatisticsByTopology {
    public static void main(String[] args) {
        Map<Boolean, StatisticsCollector.AminoAcidSummary> occurrences = MembraneConstants.PdbtmAlphaNr.getAminoAcids()
                .sequential()
                .collect(Collectors.partitioningBy(aminoAcid -> aminoAcid
                        .getParentChain()
                        .getParentProtein()
                        .getFeatureContainer()
                        .getFeature(MembraneContainer.class)
                        .isTransmembraneGroup(aminoAcid), StatisticsCollector.toAminoAcidSummary()));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("global_by_topology.tsv"),
                "header\t" + occurrences.get(true).getHeaderLine() + System.lineSeparator() +
                "tm\t" + occurrences.get(true).getOccurrenceLine() + System.lineSeparator() +
                "ntm\t" + occurrences.get(false).getOccurrenceLine());
    }
}
