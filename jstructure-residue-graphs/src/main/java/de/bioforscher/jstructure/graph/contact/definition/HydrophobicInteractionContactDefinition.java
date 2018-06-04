package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.feature.interaction.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class HydrophobicInteractionContactDefinition implements ContactDefinition {
    HydrophobicInteractionContactDefinition() {

    }

    @Override
    public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrophobicInteraction.class);
    }
}
