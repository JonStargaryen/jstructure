package studies.membrane.pdbtm;

import studies.membrane.MembraneConstants;

/**
 * Accumulate general statistics describing the data set.
 * Consider split between respectively distribution of:
 * <ul>
 *     <li>TM and nTM regions</li>
 *     <li>sequence motifs</li>
 *     <li>secondary structure elements</li>
 *     <li>interaction types</li>
 *     <li>interacting parts of amino acids (backbone, side chain, mix?)</li>
 * </ul>
 * Created by bittrich on 6/7/17.
 */
public class T01_GeneralStatistics {
    public static void main(String[] args) {
        String output = MembraneConstants.PdbtmAlphaNr.getAminoAcids()
                .sequential()
                .collect(StatisticsCollector.toAminoAcidSummary())
                .getLine();

        MembraneConstants.write(MembraneConstants.PDBTM_STATISTICS_PATH.resolve("global.tsv"), output);
    }
}
