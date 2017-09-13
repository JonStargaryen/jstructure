package de.bioforscher.jstructure.membrane.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

public class A03_RunNetCarto extends NetCartoBridge {
    private A03_RunNetCarto() {
        super(MembraneConstants.PDBTM_NR_ALPHA_DATASET_NETWORK_DIRECTORY);
    }

    public static void main(String[] args) {
        new A03_RunNetCarto();
    }
}
