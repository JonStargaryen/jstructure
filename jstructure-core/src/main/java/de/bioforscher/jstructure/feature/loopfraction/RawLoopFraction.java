package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The raw (not smoothed, binary) loop fraction values.
 * Created by bittrich on 5/17/17.
 */
public class RawLoopFraction extends FeatureContainerEntry {
    private final double rawLoopFraction;

    public RawLoopFraction(AbstractFeatureProvider featureProvider, double rawLoopFraciton) {
        super(featureProvider);
        this.rawLoopFraction = rawLoopFraciton;
    }

    public double getRawLoopFraction() {
        return rawLoopFraction;
    }
}
