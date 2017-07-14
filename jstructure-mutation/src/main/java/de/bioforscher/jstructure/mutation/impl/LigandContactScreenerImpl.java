package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.LigandContactScreener;

/**
 * Implementation of the ligand contact screener.
 * Created by bittrich on 7/13/17.
 */
public class LigandContactScreenerImpl implements LigandContactScreener {
    private final double interactionCutoff;

    public LigandContactScreenerImpl(double interactionCutoff) {
        this.interactionCutoff = interactionCutoff;
    }

    public LigandContactScreenerImpl() {
        this(DEFAULT_INTERACTION_CUTOFF);
    }

    @Override
    public double getInteractionCutoff() {
        return interactionCutoff;
    }

    @Override
    public int determineNumberOfLigandContacts(Protein protein, AminoAcid aminoAcid) {
        double[] coordinates = aminoAcid.calculate().center().getValue();
        return (int) protein
                .select()
                .ligands()
                .negationModeEnter()
                .water()
                .negationModeLeave()
                .groupDistance(coordinates, interactionCutoff)
                .asFilteredGroups()
                .count();
    }
}
