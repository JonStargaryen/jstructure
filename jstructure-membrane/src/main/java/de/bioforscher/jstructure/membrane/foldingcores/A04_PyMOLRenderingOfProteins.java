package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.membrane.MembraneConstants;
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
//            logger.info("processing {} chain {}",
//                    pdbId,
//                    chainId);

//            Structure structure = StructureParser.source(pdbId)
//                    .parse();
//            Chain chain = structure.select()
//                    .chainId(chainId)
//                    .asChain();

            List<String> residueIdentifiers;
            try (Stream<String> lines = Files.lines(path)) {
                residueIdentifiers = lines.filter(line -> !line.startsWith("#"))
                        .filter(line -> line.endsWith("EARLY"))
                        .map(line -> line.split(";")[0])
                        .collect(Collectors.toList());
            }

//            logger.info("earling folding residues are:{}{}",
//                    System.lineSeparator(),
//                    residueIdentifiers);

            Path outputPath = MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("pymol").resolve(id + ".pml");
            String pymolCommandString = residueIdentifiers.stream()
                    .map(A04_PyMOLRenderingOfProteins::handleResidueIdentifier)
                    .collect(Collectors.joining(System.lineSeparator(),
                            "reinit" + System.lineSeparator() +
                            "bg_color white" + System.lineSeparator() +
                            "unset depth_cue" + System.lineSeparator() +
                            "set ray_opaque_background, off" + System.lineSeparator() +
                            "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                            "bg_color white" + System.lineSeparator() +
                            // hide non-relevant stuff
                            "hide everything" + System.lineSeparator() +
                            "show cartoon, chain " + chainId + System.lineSeparator() +
                            // decolor everything
                            "color grey80" + System.lineSeparator() +
                            // create custom color palette
                            "set_color efr=[" + StandardFormat.format(23 / 255.0) + ", " + StandardFormat.format(111 / 255.0) + ", " + StandardFormat.format(193 / 255.0) + "]" + System.lineSeparator() +
                            "set ray_trace_mode, 3" + System.lineSeparator() +
                            "orient" + System.lineSeparator(),
                            System.lineSeparator() + "ray" + System.lineSeparator() +
                            "png " + outputPath.getParent().resolve(id + ".png") + System.lineSeparator()));

            MembraneConstants.write(outputPath, pymolCommandString);

            System.out.println("@" + outputPath);
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    pdbId,
                    e);
        }
    }

    private static String handleResidueIdentifier(String residueIdentifier) {
        return "color efr, resi " + residueIdentifier;
    }
}
