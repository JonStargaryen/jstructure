package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

public class A02_ExtractNetworks extends NetworkExtractor {
    private A02_ExtractNetworks() {
        super(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY);
    }

    public static void main(String[] args) {
        new A02_ExtractNetworks();
    }
}
