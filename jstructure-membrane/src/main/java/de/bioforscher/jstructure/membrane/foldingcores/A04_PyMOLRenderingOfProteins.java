package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A04_PyMOLRenderingOfProteins {
    private static final Logger logger = LoggerFactory.getLogger(A04_PyMOLRenderingOfProteins.class);

    public static void main(String[] args) {
        MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .forEach(A04_PyMOLRenderingOfProteins::handleProtein);
    }

    private static void handleProtein(Path path) {
        String pdbId = "?";
        try {
            String id;
            try(Stream<String> lines = Files.lines(path)) {
                id = lines.filter(line -> line.startsWith("#pdb:"))
                        .findFirst()
                        .get()
                        .split(": ")[1];
            }
            pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            logger.info("processing {} chain {}",
                    pdbId,
                    chainId);

            Structure structure = StructureParser.source(pdbId)
                    .parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            List<String> residueIdentifiers;
            try (Stream<String> lines = Files.lines(path)) {
                residueIdentifiers = lines.filter(line -> !line.startsWith("#"))
                        .filter(line -> line.endsWith("EARLY"))
                        .map(line -> line.split(";")[0])
                        .collect(Collectors.toList());
            }

            logger.info("earling folding residues are:{}{}",
                    System.lineSeparator(),
                    residueIdentifiers);
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    pdbId,
                    e);
        }
    }
}
