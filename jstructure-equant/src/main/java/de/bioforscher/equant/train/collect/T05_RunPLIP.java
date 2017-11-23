package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.jstructure.feature.interactions.PLIPRestServiceQuery;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Query the PLIPRestService for interaction annotations of predicted models.
 */
public class T05_RunPLIP {
    private static final Logger logger = LoggerFactory.getLogger(T05_RunPLIP.class);
    private static final Path DIRECTORY = EquantConstants.CASP9_DIRECTORY.resolve("predictions");
    private static final Path OUTPUT_DIRECTORY = EquantConstants.CASP9_DIRECTORY.resolve("plip");

    public static void main(String[] args) {
        EquantConstants.list(DIRECTORY)
                .forEach(T05_RunPLIP::handleDirectory);
    }

    private static void handleDirectory(Path inDirectory) {
        Path outDirectory = OUTPUT_DIRECTORY.resolve(inDirectory.toFile().getName());
        EquantConstants.list(inDirectory)
                .forEach(path -> {
                    try {
                        logger.info("processing {}",
                                path.toFile().getName());

                        String name = path.toFile().getName();
                        Path outputPath = outDirectory.resolve(name + ".plip");
                        if(Files.exists(outputPath)) {
                            if(EquantConstants.lines(outputPath)
                                    .anyMatch(line -> line.contains("has_interactions=\"True\""))) {
                                logger.info("skipping already processed, sane file");
                                return;
                            }
                        }

                        Structure structure = StructureParser.source(path).parse();
                        Chain chain = structure.chains().findFirst().get();
                        Document document = PLIPRestServiceQuery.calculateIntraChainDocument(chain);
                        EquantConstants.write(outputPath, document.html());
                    } catch (Exception e) {
                        logger.warn("failed to process {}",
                                path.toFile().getName(),
                                e);
                    }
                });
    }
}
