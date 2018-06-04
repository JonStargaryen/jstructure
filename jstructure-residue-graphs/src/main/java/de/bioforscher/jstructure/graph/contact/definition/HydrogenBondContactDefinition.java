package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.feature.interaction.HydrogenBond;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class HydrogenBondContactDefinition implements ContactDefinition {
    HydrogenBondContactDefinition() {

    }

    @Override
    public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrogenBond.class);
    }
}
