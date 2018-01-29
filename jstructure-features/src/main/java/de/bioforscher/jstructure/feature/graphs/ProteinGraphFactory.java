package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.feature.interactions.HydrogenBond;
import de.bioforscher.jstructure.feature.interactions.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;

public class ProteinGraphFactory {
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

    public static ProteinGraph createProteinGraph(Chain chain, InteractionScheme interactionScheme) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Pair<AminoAcid, AminoAcid>> contacts = new ArrayList<>();

        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                    continue;
                }

                if(interactionScheme.areInContact(aminoAcid1, aminoAcid2)) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                }
            }
        }

        return new ProteinGraph(aminoAcids, contacts);
    }
}
