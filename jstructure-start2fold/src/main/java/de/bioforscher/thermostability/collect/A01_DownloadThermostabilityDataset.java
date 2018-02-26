package de.bioforscher.thermostability.collect;

import de.bioforscher.start2fold.Start2FoldConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class A01_DownloadThermostabilityDataset {
    private static final Logger logger = LoggerFactory.getLogger(A01_DownloadThermostabilityDataset.class);
    private static final Path PDB_DIRECTORY = Start2FoldConstants.THERMOSTABILITY_DIRECTORY.resolve("pdb");

    public static void main(String[] args) throws IOException {
        List<String> lines = Files.readAllLines(Start2FoldConstants.THERMOSTABILITY_DIRECTORY.resolve("sadeghi.list"));

        // download PDB structures
        for (String line : lines) {
            String[] split = line.split(";");
            String thermophilePdbId = split[0].toLowerCase();
            String mesophilePdbId = split[2].toLowerCase();

            try {
                logger.info("downloading thermophile {}", thermophilePdbId);
                Start2FoldConstants.move(new URL("https://files.rcsb.org/download/" + thermophilePdbId + ".pdb"),
                        PDB_DIRECTORY.resolve(thermophilePdbId + ".pdb"));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            try {
                logger.info("downloading mesophile {}", mesophilePdbId);
                Start2FoldConstants.move(new URL("https://files.rcsb.org/download/" + mesophilePdbId + ".pdb"),
                        PDB_DIRECTORY.resolve(mesophilePdbId + ".pdb"));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }
}
