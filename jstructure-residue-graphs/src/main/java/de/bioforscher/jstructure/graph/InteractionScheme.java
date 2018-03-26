package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.feature.interaction.HydrogenBond;
import de.bioforscher.jstructure.feature.interaction.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.function.BiPredicate;

public enum InteractionScheme {
    CALPHA8((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
            .distanceFast(aminoAcid2.getCa()) < 8 * 8),
    SALENTIN2015((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
            .areInContact(aminoAcid1, aminoAcid2)),
    SALENTIN2015_HYDROGEN_BONDS(((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
            .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrogenBond.class))),
    SALENTIN2015_HYDROPHOBIC_INTERACTION(((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
            .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrophobicInteraction.class)));

    private final BiPredicate<AminoAcid, AminoAcid> criterion;

    InteractionScheme(BiPredicate<AminoAcid, AminoAcid> criterion) {
        this.criterion = criterion;
    }

    public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return criterion.test(aminoAcid1, aminoAcid2);
    }
}
