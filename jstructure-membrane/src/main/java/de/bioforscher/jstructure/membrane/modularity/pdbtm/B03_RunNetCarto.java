package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

public class B03_RunNetCarto extends NetCartoBridge {
    private B03_RunNetCarto() {
        super(MembraneConstants.PDBTM_NR_BETA_DATASET_NETWORK_DIRECTORY);
    }

    public static void main(String[] args) {
        new B03_RunNetCarto();
    }
}
