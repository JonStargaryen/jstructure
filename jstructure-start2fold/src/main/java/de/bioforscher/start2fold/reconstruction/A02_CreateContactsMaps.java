package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.start2fold.Start2FoldConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.NoSuchElementException;

public class A02_CreateContactsMaps {
    private static final Logger logger = LoggerFactory.getLogger(A02_CreateContactsMaps.class);
    private static final Path DIRECTORY = Start2FoldConstants.START2FOLD_DIRECTORY.resolve("maps");

    public static void main(String[] args) {
        Start2FoldConstants.lines(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("ids.list"))
                .forEach(A02_CreateContactsMaps::handlePdbId);
    }

    private static void handlePdbId(String pdbId) {
        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();
        logger.info("processing {} with {} chains",
                pdbId,
                structure.chains().count());

        Path start2foldPath = Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("start2fold"))
                .filter(path -> {
                    try {
                        return Files.lines(path)
                                .anyMatch(line -> line.contains(pdbId));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                })
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("did not find file describing " + pdbId));

        new ReconstructionContactMapCreator(chain, DIRECTORY, start2foldPath);
    }
}
