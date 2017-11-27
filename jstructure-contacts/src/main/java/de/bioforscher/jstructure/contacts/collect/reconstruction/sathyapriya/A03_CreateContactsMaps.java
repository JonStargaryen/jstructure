package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.ContactMapCreator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;

public class A03_CreateContactsMaps {
    private static final Logger logger = LoggerFactory.getLogger(A03_CreateContactsMaps.class);
    private static final Path DIRECTORY = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("maps");

    public static void main(String[] args) {
        ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .forEach(A03_CreateContactsMaps::handlePdbId);
    }

    private static void handlePdbId(String pdbId) {
        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();
        logger.info("processing {} with {} chains",
                pdbId,
                structure.chains().count());

        new ContactMapCreator(chain, DIRECTORY);
    }
}
