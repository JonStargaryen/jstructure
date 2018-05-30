package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.feature.interaction.HydrogenBond;
import de.bioforscher.jstructure.feature.interaction.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.SelectionException;

import java.util.function.BiPredicate;

public enum InteractionScheme {
    //TODO changeable contact definition
    CALPHA10((aminoAcid1, aminoAcid2) -> {
        try {
            return aminoAcid1.getCa().calculate()
                    .distanceFast(aminoAcid2.getCa()) < 100;
        } catch (NullPointerException | SelectionException e) {
            return aminoAcid1.calculate().centroid()
                    .distanceFast(aminoAcid2.calculate().centroid()) < 100;
        }
    }),
    CBETA9((aminoAcid1, aminoAcid2) -> {
        try {
            return aminoAcid1.select()
                    .betaCarbonAtoms()
                    .asAtom().calculate()
                    .distanceFast(aminoAcid2.select()
                            .betaCarbonAtoms()
                            .asAtom()) < 9 * 9;
        } catch (NullPointerException | SelectionException e) {
            return aminoAcid1.calculate().centroid()
                    .distanceFast(aminoAcid2.calculate().centroid()) < 9 * 9;
        }
    }),
    CALPHA8((aminoAcid1, aminoAcid2) -> {
        try {
            return aminoAcid1.select()
                    .betaCarbonAtoms()
                    .asAtom().calculate()
                    .distanceFast(aminoAcid2.select()
                            .betaCarbonAtoms()
                            .asAtom()) < 9 * 9;
        } catch (NullPointerException | SelectionException e) {
            return aminoAcid1.calculate().centroid()
                    .distanceFast(aminoAcid2.calculate().centroid()) < 9 * 9;
        }
    }),
//    CALPHA8((aminoAcid1, aminoAcid2) -> {
//        try {
//            return aminoAcid1.getCa().calculate()
//                    .distanceFast(aminoAcid2.getCa()) < 8 * 8;
//        } catch (NullPointerException | SelectionException e) {
//            return aminoAcid1.calculate().centroid()
//                    .distanceFast(aminoAcid2.calculate().centroid()) < 8 * 8;
//        }
//    }),
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
