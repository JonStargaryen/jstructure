package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

public class ContactMapScorer {
    private static final Logger logger = LoggerFactory.getLogger(ContactMapScorer.class);

    public static ContactMapScoring score(Structure reference, Structure model) {
        logger.info("creating contact map representation of {} and {}",
                reference.getProteinIdentifier().getFullName(),
                model.getProteinIdentifier().getFullName());

        StructureRenumberer.renumberStructure(reference, model);

        List<AminoAcid> aminoAcids = reference.aminoAcids()
                .collect(Collectors.toList());

        ContactMap plipReferenceContactMap = ContactMap.createPlipContactMap(reference);
        ContactMap plipModelContactMap = ContactMap.createPlipContactMap(model);

        ContactMap conventionalReferenceContactMap = ContactMap.createBetaCarbon8AngstromContactMap(reference);
        ContactMap conventionalModelContactMap = ContactMap.createBetaCarbon8AngstromContactMap(model);

        return new ContactMapScoring(aminoAcids,
                conventionalReferenceContactMap,
                conventionalModelContactMap,
                plipReferenceContactMap,
                plipModelContactMap);
    }
}
