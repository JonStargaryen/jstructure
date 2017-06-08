package studies.membrane.pdbtm;

import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import studies.membrane.MembraneConstants;

import java.util.Map;
import java.util.stream.Collectors;

/**
 * Statistics now with a split by topology.
 * Created by bittrich on 6/7/17.
 */
public class T02_StatisticsByTopology {
    public static void main(String[] args) {
        Map<Boolean, StatisticsCollector.OccurrenceSummary> occurrences = MembraneConstants.getAminoAcidsOfPdbtmAlphaNrList()
                .sequential()
                .collect(Collectors.partitioningBy(aminoAcid -> aminoAcid
                        .getParentChain()
                        .getParentProtein()
                        .getFeatureContainer()
                        .getFeature(MembraneContainer.class)
                        .isTransmembraneGroup(aminoAcid), StatisticsCollector.toOccurrenceSummary()));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("global_by_topology.tsv"),
                "header\t" + StatisticsCollector.OccurrenceSummary.getHeaderLine() + System.lineSeparator() +
                "tm\t" + occurrences.get(true).getOccurrenceLine() + System.lineSeparator() +
                "ntm\t" + occurrences.get(false).getOccurrenceLine());
    }
}
