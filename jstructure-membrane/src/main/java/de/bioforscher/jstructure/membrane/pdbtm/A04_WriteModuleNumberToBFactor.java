package de.bioforscher.jstructure.membrane.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;

public class A04_WriteModuleNumberToBFactor extends ModuleToBFactorWriter {
    private A04_WriteModuleNumberToBFactor(Path datasetPath) {
        super(datasetPath);
    }

    public static void main(String[] args) {
        new A04_WriteModuleNumberToBFactor(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY);
    }
}