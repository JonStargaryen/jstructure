package de.bioforscher.jstructure.membrane.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;

public class B01_DownloadDataset extends PdbtmDatasetComposer {
    private static final String FETCH_URL = "http://pdbtm.enzim.hu/data/pdbtm_beta_nr.list";

    private B01_DownloadDataset(String fetchUrl, Path outputPath) {
        super(fetchUrl, outputPath);
    }

    public static void main(String[] args) {
        new B01_DownloadDataset(FETCH_URL, MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY);
    }
}
