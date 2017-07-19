package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * Reports the conservation according to some criteria for a given amino acid.
 * Created by bittrich on 7/13/17.
 */
public class ConservationProfile extends FeatureContainerEntry {
    private final double sequentialConservation;
    private final double energeticalConservation;

    public ConservationProfile(double sequence, double energy) {
        super(null);
        this.sequentialConservation = sequence;
        this.energeticalConservation = energy;
    }

    public double getSequentialConservation() {
        return sequentialConservation;
    }

    public double getEnergeticalConservation() {
        return energeticalConservation;
    }
}
