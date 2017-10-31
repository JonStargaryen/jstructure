package de.bioforscher.jstructure.membrane.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.DatasetComposer;

import java.nio.file.Path;

public class A01_DownloadDataset extends DatasetComposer {
    private static final String FETCH_URL = "http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list";

    private A01_DownloadDataset(String fetchUrl, Path outputPath) {
        super(fetchUrl, outputPath);
    }

    public static void main(String[] args) {
        new A01_DownloadDataset(FETCH_URL, MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY);
    }
}
