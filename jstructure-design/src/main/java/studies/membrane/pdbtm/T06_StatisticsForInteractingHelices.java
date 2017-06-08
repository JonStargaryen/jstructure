package studies.membrane.pdbtm;

import studies.membrane.MembraneConstants;

/**
 * Extract interactions between amino acids in different helices.
 * Created by bittrich on 6/8/17.
 */
public class T06_StatisticsForInteractingHelices {
    public static void main(String[] args) {
        StatisticsCollector.InteractingAminoAcidSummary summary = MembraneConstants.PdbtmAlphaNr
                .getInteractionsBetweenTransmembraneHelices()
                .sequential()
                .collect(StatisticsCollector.toInteractingAminoAcidSummary());

        String output = summary.getLine();

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("helixInteractions.tsv"), output);
    }
}
