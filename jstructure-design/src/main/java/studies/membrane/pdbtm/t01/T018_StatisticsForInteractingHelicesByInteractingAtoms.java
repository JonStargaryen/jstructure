package studies.membrane.pdbtm.t01;

import studies.membrane.MembraneConstants;
import studies.membrane.pdbtm.StatisticsCollector;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Extract interactions between amino acids in different helices, grouped by interacting atoms.
 * Created by bittrich on 6/8/17.
 */
public class T018_StatisticsForInteractingHelicesByInteractingAtoms {
    public static void main(String[] args) {
        Map<String, StatisticsCollector.InteractingAminoAcidSummary> statistics = new HashMap<>();
        statistics.put("backbone", new StatisticsCollector.InteractingAminoAcidSummary());
        statistics.put("sideChain", new StatisticsCollector.InteractingAminoAcidSummary());
        statistics.put("mixed", new StatisticsCollector.InteractingAminoAcidSummary());

        MembraneConstants.PdbtmAlphaNr.getInteractionsBetweenTransmembraneHelices()
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
                    StatisticsCollector.InteractingAminoAcidSummary aminoAcidSummary = statistics.get(atoms);
                    aminoAcidSummary.accept(plipInteraction);
                });

        String output = statistics.entrySet().stream()
                .map(entry -> entry.getKey() + "\t" + entry.getValue().getOccurrenceLine())
                .collect(Collectors.joining(System.lineSeparator(),
                        "atoms\t" + statistics.values().iterator().next().getHeaderLine() + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("interactingHelices_by_interactingAtoms.tsv"), output);
    }
}
