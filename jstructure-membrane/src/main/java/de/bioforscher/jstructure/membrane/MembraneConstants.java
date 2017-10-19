package de.bioforscher.jstructure.membrane;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

/**
 * Shared class of project constants and convenience functions used in the membrane modules.
 */
public class MembraneConstants {
    private static final Logger logger = LoggerFactory.getLogger(MembraneConstants.class);
    public static final Path GIT_DIRECTORY = Paths.get(System.getProperty("user.home")).resolve("git");
    public static final Path PHD_DIRECTORY = GIT_DIRECTORY.resolve("phd_sb_repo");
    public static final Path DATA_DIRECTORY = PHD_DIRECTORY.resolve("data");
    public static final Path DATASETS_DIRECTORY = DATA_DIRECTORY.resolve("datasets");
    // test performance and optimal setup of module annotation
    public static final Path MODULARITY_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("modularity");
    public static final Path PDBTM_NR_ALPHA_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("pdbtm_nr_alpha");
    public static final Path PDBTM_NR_ALPHA_DATASET_PDB_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("pdb");
    public static final Path PDBTM_NR_ALPHA_DATASET_OPM_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("opm");
    public static final Path PDBTM_NR_ALPHA_DATASET_PLIP_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("plip");
    public static final Path PDBTM_NR_ALPHA_DATASET_NETWORK_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("network");
    public static final Path PDBTM_NR_BETA_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("pdbtm_nr_beta");
    public static final Path PDBTM_NR_BETA_DATASET_PDB_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("pdb");
    public static final Path PDBTM_NR_BETA_DATASET_OPM_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("opm");
    public static final Path PDBTM_NR_BETA_DATASET_PLIP_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("plip");
    public static final Path PDBTM_NR_BETA_DATASET_NETWORK_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("network");
    public static final Path FOLDING_CORES_DIRECTORY = DATASETS_DIRECTORY.resolve("foldingcores");

    public static Stream<String> lines(Path path) {
        try {
            return Files.lines(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<String> lines(String url) {
        try {
            return new BufferedReader(new InputStreamReader(new URL(url).openStream())).lines();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<Path> list(Path path) {
        try {
            return Files.list(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, String content) {
        ensureDirectoriesExist(path);
        try {
            Files.write(path, content.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void ensureDirectoriesExist(Path path) {
        if(Files.isDirectory(path)) {
            throw new IllegalArgumentException(path.toFile().getAbsolutePath() + " is no regular file - cannot process");
        }

        Path parentDirectory = path.getParent();
        if(!Files.exists(parentDirectory)) {
            logger.info("creating directory '{}'", parentDirectory.toFile().getAbsolutePath());
            try {
                Files.createDirectories(parentDirectory);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }
}
