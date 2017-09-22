package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;

public class B04_WriteModuleNumberToBFactor extends ModuleToBFactorWriter {
    private B04_WriteModuleNumberToBFactor(Path datasetPath) {
        super(datasetPath);
    }

    public static void main(String[] args) {
        new B04_WriteModuleNumberToBFactor(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY);
    }
}
