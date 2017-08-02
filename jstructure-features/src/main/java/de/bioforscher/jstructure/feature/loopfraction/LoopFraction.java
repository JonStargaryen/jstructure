package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

/**
 * The smoothed loop fraction for groups.
 * Created by bittrich on 5/17/17.
 */
public class LoopFraction extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double loopFraction;

    LoopFraction(FeatureProvider featureProvider, double loopFraction) {
        super(featureProvider);
        this.loopFraction = loopFraction;
    }

    public double getLoopFraction() {
        return loopFraction;
    }

    @Override
    public Double getValue() {
        return getLoopFraction();
    }
}
