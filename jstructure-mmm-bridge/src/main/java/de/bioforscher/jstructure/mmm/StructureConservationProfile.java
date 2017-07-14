package de.bioforscher.jstructure.mmm;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

/**
 * Reports the structural conservation according to the number of mmm hits for a given amino acid. Wraps the structural
 * conservation score as feature.
 * Created by bittrich on 7/13/17.
 */
public class StructureConservationProfile extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double structuralConservationScore;

    public StructureConservationProfile(double value) {
        super(null);
        this.structuralConservationScore = value;
    }

    @Override
    public Double getValue() {
        return structuralConservationScore;
    }

    public double getStructuralConservationScore() {
        return structuralConservationScore;
    }
}
