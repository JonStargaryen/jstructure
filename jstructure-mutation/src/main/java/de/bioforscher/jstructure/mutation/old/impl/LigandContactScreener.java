package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Screens for ligands around a given amino acid.
 * Created by bittrich on 7/12/17.
 */
@Deprecated
public class LigandContactScreener {
    private static final double INTERACTION_CUTOFF = 4.0;

    public static double determineNumberOfLigandContacts(Group group) {
        return group.getParentChain()
                .select()
                .ligands()
                .negationModeEnter()
                .water()
                .negationModeLeave()
                .groupDistance(group.calculate().centroid().getValue(), INTERACTION_CUTOFF)
                .asFilteredGroups()
                .count();
    }
}
