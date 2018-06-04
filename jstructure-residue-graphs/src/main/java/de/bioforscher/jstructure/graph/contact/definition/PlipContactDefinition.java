package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class PlipContactDefinition implements ContactDefinition {
    PlipContactDefinition() {

    }

    @Override
    public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContact(aminoAcid1, aminoAcid2);
    }
}
