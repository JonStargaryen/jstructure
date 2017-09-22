package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

public class B02_ExtractNetworks extends NetworkExtractor {
    private B02_ExtractNetworks() {
        super(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY);
    }

    public static void main(String[] args) {
        new B02_ExtractNetworks();
    }
}
