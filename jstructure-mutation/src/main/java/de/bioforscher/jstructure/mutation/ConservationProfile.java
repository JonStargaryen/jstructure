package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Reports the conservation according to some criteria for a given amino acid.
 * Created by bittrich on 7/13/17.
 */
public class ConservationProfile extends FeatureContainerEntry {
    private final double sequentialConservation;
    private final double structuralConservation;
    private final double energeticalConservation;

    public ConservationProfile(double sequence, double structure, double energy) {
        super(null);
        this.sequentialConservation = sequence;
        this.structuralConservation = structure;
        this.energeticalConservation = energy;
    }

    public double getSequentialConservation() {
        return sequentialConservation;
    }

    public double getStructuralConservation() {
        return structuralConservation;
    }

    public double getEnergeticalConservation() {
        return energeticalConservation;
    }
}
