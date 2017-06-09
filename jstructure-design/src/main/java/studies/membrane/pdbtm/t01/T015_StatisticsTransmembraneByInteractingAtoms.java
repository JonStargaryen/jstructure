package studies.membrane.pdbtm.t01;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import studies.membrane.MembraneConstants;
import studies.membrane.pdbtm.StatisticsCollector;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Statistics for transmembrane regions, grouped by interacting atoms.
 * Created by bittrich on 6/8/17.
 */
public class T015_StatisticsTransmembraneByInteractingAtoms {
    public static void main(String[] args) {
        Map<String, StatisticsCollector.AminoAcidSummary> statistics = new HashMap<>();
        statistics.put("backbone", new StatisticsCollector.AminoAcidSummary());
        statistics.put("sideChain", new StatisticsCollector.AminoAcidSummary());
        statistics.put("mixed", new StatisticsCollector.AminoAcidSummary());

        MembraneConstants.PdbtmAlphaNr.getInteractionsTransmembrane()
                .sequential()
                .forEach(plipInteraction -> {
                    String atoms;
                    if(plipInteraction.isBackboneInteraction()) {
                        atoms = "backbone";
                    } else if(plipInteraction.isSideChainInteraction()) {
                        atoms = "sideChain";
                    } else {
                        atoms = "mixed";
                    }
                    StatisticsCollector.AminoAcidSummary aminoAcidSummary = statistics.get(atoms);
                    aminoAcidSummary.accept((AminoAcid) plipInteraction.getPartner1());
                    aminoAcidSummary.accept((AminoAcid) plipInteraction.getPartner2());
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "atoms\t" + statistics.values().iterator().next().getHeaderLine() + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("transmembrane_by_interactingAtoms.tsv"), output);
    }
}
