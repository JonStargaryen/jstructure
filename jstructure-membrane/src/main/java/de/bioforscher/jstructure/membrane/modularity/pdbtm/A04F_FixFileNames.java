package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 * Fixes wrong file names (which were introduced by shitty batch-scripting).
 */
public class A04F_FixFileNames {
    private static final Logger logger = LoggerFactory.getLogger(A04F_FixFileNames.class);

    public static void main(String[] args) {
        MembraneConstants.list(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("network"))
                .filter(path -> path.toFile().getName().contains("_plip.dat."))
//        MembraneConstants.list(MembraneConstants.MODULARITY_DATASET_DIRECTORY.resolve("network"))
//                .filter(path -> path.toFile().getName().contains(".dat."))
                .forEach(path -> {
                    try {
                        String filename = path.toFile().getName();
                        String[] split = filename.split("\\.");
                        String correctFilename = split[0] + "." + split[2] + "." + split[3];
                        logger.info("renaming {} to {}",
                                filename,
                                correctFilename);
                        Files.move(path, path.getParent().resolve(correctFilename), StandardCopyOption.REPLACE_EXISTING);
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}
