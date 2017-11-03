package de.bioforscher.jstructure.membrane.collect.kinks;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Files;
import java.nio.file.Paths;

public class A03_MoveKinkFinderResults {
    public static void main(String[] args) {
        MembraneConstants.list(Paths.get("/home/bittrich/Downloads/KF_err_lin/"))
                .filter(Files::isDirectory)
                .filter(path -> path.toFile().getName().contains("results_"))
                .map(path -> path.resolve("kinks.csv"))
                .forEach(path -> {
                    String pdbId = path.getParent().toFile().getName().split("_")[1];
                    System.out.println(pdbId);
                    MembraneConstants.move(path, MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY
                            .resolve("kinks")
                            .resolve(pdbId + ".kinks"));
                });
    }
}
