package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

/**
 * The raw (not smoothed, binary) loop fraction values.
 * Created by bittrich on 5/17/17.
 */
@DefaultFeatureProvider(LoopFractionCalculator.class)
public class RawLoopFraction extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double rawLoopFraction;

    RawLoopFraction(FeatureProvider featureProvider, double rawLoopFraciton) {
        super(featureProvider);
        this.rawLoopFraction = rawLoopFraciton;
    }

    public double getRawLoopFraction() {
        return rawLoopFraction;
    }

    @Override
    public Double getValue() {
        return getRawLoopFraction();
    }
}
