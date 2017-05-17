package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The smoothed loop fraction for groups.
 * Created by bittrich on 5/17/17.
 */
public class LoopFraction extends FeatureContainerEntry {
    private final double loopFraction;

    public LoopFraction(AbstractFeatureProvider featureProvider, double loopFraction) {
        super(featureProvider);
        this.loopFraction = loopFraction;
    }

    public double getLoopFraction() {
        return loopFraction;
    }
}
