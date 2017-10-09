package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.modularity.pdbtm.ModuleToBFactorWriter;
import de.bioforscher.jstructure.membrane.modularity.pdbtm.NetworkExtractor;
import de.bioforscher.jstructure.membrane.modularity.pdbtm.PdbtmDatasetComposer;

import java.nio.file.Path;

/**
 * Compare NetCarto modules and parameters with experimental data to ensure a valid setup.
 * Especially wants to check:
 * <ul>
 *     <li>the best contact definition: PLIP, 0.6 A (Krishnan, 2007), 1.0 A (Khan, 2015)</li>
 *     <li>whether consecutive amino acids should be linked explicitly</li>
 * </ul>
 *
 * Available structures:
 * <ul>
 *     <li>Hleap, 2013: Isoamylase (PDB: 1bf2_A)</li>
 *     <li>Englander, 2014 -> Hu, 2013: Ribonuclease H (PDB: 1f21_A)</li>
 *     <li>Englander, 2014 -> Roder, 1988 -> Takano, 1981, Bai, 1995: Cyt c (PDB: 3cyt_I ?)</li>
 *     <li>Hleap, 2013: Niemann-Pick disease, type C1 protein (PDB: 3gkh_A)</li>
 * </ul>
 */
public class NetCartoValidation {
    public static void main(String[] args) {
        Path dataset = MembraneConstants.MODULARITY_DATASET_DIRECTORY;

        // create dataset files from definition file
        new PdbtmDatasetComposer(dataset);

        new NetworkExtractor(dataset);

        new ModuleToBFactorWriter(dataset);
    }
}
