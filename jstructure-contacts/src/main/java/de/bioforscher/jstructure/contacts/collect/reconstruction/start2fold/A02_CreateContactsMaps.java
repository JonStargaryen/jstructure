package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.ReconstructionContactMapCreator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.NoSuchElementException;

public class A02_CreateContactsMaps {
    private static final Logger logger = LoggerFactory.getLogger(A02_CreateContactsMaps.class);
    private static final Path DIRECTORY = ContactsConstants.START2FOLD_DIRECTORY.resolve("maps");

    public static void main(String[] args) {
        ContactsConstants.lines(ContactsConstants.START2FOLD_DIRECTORY.resolve("ids.list"))
                .forEach(A02_CreateContactsMaps::handlePdbId);
    }

    private static void handlePdbId(String pdbId) {
        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();
        logger.info("processing {} with {} chains",
                pdbId,
                structure.chains().count());

        Path start2foldPath = ContactsConstants.list(ContactsConstants.START2FOLD_DIRECTORY.resolve("start2fold"))
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
