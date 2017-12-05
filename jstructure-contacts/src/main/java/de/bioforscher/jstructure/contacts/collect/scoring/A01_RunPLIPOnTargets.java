package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.PLIPRestServiceQuery;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;

public class A01_RunPLIPOnTargets {
    private static final Logger logger = LoggerFactory.getLogger(A01_RunPLIPOnTargets.class);
    private static final Path BASE_DIRECTORY = ContactsConstants.GIT_DIRECTORY.resolve("phd_sb_repo")
            .resolve("data")
            .resolve("equant")
            .resolve("casp")
            .resolve("CASP9");
    private static final Path TARGET_DIRECTORY = BASE_DIRECTORY.resolve("targets");
    private static final Path OUTPUT_DIRECTORY = BASE_DIRECTORY.resolve("plip");

    public static void main(String[] args) {
        ContactsConstants.list(TARGET_DIRECTORY)
                .forEach(A01_RunPLIPOnTargets::handleTarget);
    }

    private static void handleTarget(Path targetPath) {
        Path outDirectory = OUTPUT_DIRECTORY.resolve(targetPath.toFile().getName().split("\\.")[0]);

        logger.info("processing {}",
                targetPath.toFile().getName());

        try {
            String name = targetPath.toFile().getName().split("\\.")[0];
            Path outputPath = outDirectory.resolve(name + ".plip");
            if (Files.exists(outputPath)) {
                if (ContactsConstants.lines(outputPath)
                        .anyMatch(line -> line.contains("has_interactions=\"True\""))) {
                    logger.info("skipping already processed, sane file");
                    return;
                }
            }

            Structure structure = StructureParser.source(targetPath).parse();
            Chain chain = structure.chains().findFirst().get();
            Document document = PLIPRestServiceQuery.calculateIntraChainDocument(chain);
            ContactsConstants.write(outputPath, document.html());
        } catch (Exception e) {
            logger.warn("computation failed",
                    e);
        }
    }
}
