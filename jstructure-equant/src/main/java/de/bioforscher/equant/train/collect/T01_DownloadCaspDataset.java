package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Download the CASP 9 dataset used for training.
 * Encompasses:
 * <ul>
 *     <li>target structures</li>
 *     <li>up to 30 predictions which ought to have a GDT >20 to ensure reasonable predictions</li>
 *     <li>the result files describing the scoring for these predictions</li>
 * </ul>
 */
public class T01_DownloadCaspDataset {
    private static final Logger logger = LoggerFactory.getLogger(T01_DownloadCaspDataset.class);
    private static final double THRESHOLD = 20.00;
    private static final int NUMBER_OF_SELECTED_TARGETS = 30;

    enum CASPRun {
        CASP9("http://www.predictioncenter.org/download_area/CASP9/",
                "http://www.predictioncenter.org/download_area/CASP9/results_LGA_sda/",
                515,
                643,
                "http://www.predictioncenter.org/download_area/CASP9/predictions/",
                "http://www.predictioncenter.org/download_area/CASP9/targets/casp9.targ_unsplit.tgz");

        private final String id;
        private final String sda;
        private final String predictions;
        private final String targets;
        private final int start;
        private final int end;

        CASPRun(String id, String sda, int start, int end, String predictions, String targets) {
            this.id = id;
            this.sda = sda;
            this.start = start;
            this.end = end;
            this.predictions = predictions;
            this.targets = targets;
        }
    }

    public static void main(String[] args) {
        processCASPRun(CASPRun.CASP9);
    }

    private static void processCASPRun(CASPRun run) {
        logger.info("handling {}", run.name());

        Path baseDir = EquantConstants.CASP_DATA_DIRECTORY.resolve(run.name());
        logger.info("writing files to {}",
                baseDir);

        Path targetDir = baseDir.resolve("targets");
        Path predictionDir = baseDir.resolve("predictions");
        Path resultsDir = baseDir.resolve("results");

        try {
            // handle target list
            Path targetArchive = targetDir.resolve("targ_unsplit.tgz");
            EquantConstants.move(new URL(run.targets), targetArchive);
            EquantConstants.extractTarGz(targetArchive);
            EquantConstants.delete(targetArchive);

            for (int i = run.start; i <= run.end; i++) {
                String targetName = "T0" + i;
                logger.info("target {}",
                        targetName);
                String sdaURL = targetName + ".SUMMARY.lga_sda" + (run.name().equals("CASP10") || run.name().equals("CASP11") ? "" : "_4") + ".txt";
                logger.info("fetching SDA file {}",
                        sdaURL);
                try {
                    List<String> models = new BufferedReader(new InputStreamReader(new URL(run.sda + sdaURL).openStream()))
                            .lines()
                            .filter(line -> !line.startsWith("NAME") && !line.contains("AL"))
                            .filter(line -> Double.valueOf(line.substring(64, 74)) > THRESHOLD)
                            .map(line -> line.split(":")[0])
                            .map(line -> line.split("\\.")[0])
                            .collect(Collectors.toList());
                    logger.info("{} models are suitable",
                            models.size());

                    Collections.shuffle(models);

                    List<String> selectedModels = models.subList(0, models.size() > NUMBER_OF_SELECTED_TARGETS ? NUMBER_OF_SELECTED_TARGETS : models.size());
                    if (selectedModels.isEmpty()) {
                        String targetToDelete = targetName + ".pdb";
                        logger.info("selected nothing due to GDT cutoff - deleting target {}",
                                targetToDelete);

                        // delete target file if not applicable
                        EquantConstants.delete(targetDir.resolve(targetToDelete));
                    } else {
                        logger.info("selected {}: {}",
                                NUMBER_OF_SELECTED_TARGETS,
                                selectedModels);

                        // download and keep selected predictions
                        String predictionArchiveName = targetName + ".tar.gz";
                        Path predictionArchive = predictionDir.resolve(predictionArchiveName);
                        EquantConstants.move(new URL(run.predictions + predictionArchiveName), predictionArchive);
                        EquantConstants.extractTarGz(predictionArchive);
                        EquantConstants.delete(predictionArchive);
                        EquantConstants.list(predictionDir.resolve(targetName))
                                .filter(model -> !selectedModels.contains(model.toFile().getName().split("\\.")[0]))
                                .forEach(EquantConstants::delete);

                        // download and keep selected results
                        String resultsArchiveName = targetName + ".tgz";
                        Path resultsArchive = resultsDir.resolve(resultsArchiveName);
                        EquantConstants.move(new URL(run.sda + resultsArchiveName), resultsArchive);
                        EquantConstants.extractTarGz(resultsArchive);
                        EquantConstants.delete(resultsArchive);
                        EquantConstants.list(resultsDir.resolve(targetName))
                                .filter(model -> !selectedModels.contains(model.toFile().getName().split("\\.")[0]))
                                .forEach(EquantConstants::delete);
                    }
                } catch (IOException e) {
                    logger.warn("failed to download {} {}.tar.gz",
                            run.predictions,
                            targetName);
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
