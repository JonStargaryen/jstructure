package de.bioforscher.jstructure.graph.contact.definition;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.SelectionException;

public class AlphaCarbonContactDefinition extends DistanceBasedContactDefinition {
    AlphaCarbonContactDefinition(double distance) {
        super(distance);
    }

    @Override
    public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        try {
            return aminoAcid1.getCa().calculate()
                    .distanceFast(aminoAcid2.getCa()) < squaredDistance;
        } catch (NullPointerException | SelectionException e) {
            return aminoAcid1.calculate().centroid()
                    .distanceFast(aminoAcid2.calculate().centroid()) < squaredDistance;
        }
    }

    @Override
    public String getConfoldRRType() {
        return "ca";
    }
}
