package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.model.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ContactMapScorer {
    private static final Logger logger = LoggerFactory.getLogger(ContactMapScorer.class);

    public static ContactMapScoring score(Structure reference, Structure model) {
        logger.info("creating contact map representation of {} and {}",
                reference.getProteinIdentifier().getFullName(),
                model.getProteinIdentifier().getFullName());

        StructureRenumberer.renumberStructure(reference, model);

        ContactMap plipReferenceContactMap = ContactMap.createPlipContactMap(reference);
        ContactMap plipModelContactMap = ContactMap.createPlipContactMap(model);

        ContactMap conventionalReferenceContactMap = ContactMap.createBetaCarbon8AngstromContactMap(reference);
        ContactMap conventionalModelContactMap = ContactMap.createBetaCarbon8AngstromContactMap(model);



        return null;
    }
}
