package de.bioforscher.equant.train;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Structure;

/**
 * Smooths the features used in eQuant with an arbitrary window size. Creates a new instance of FeatureContainerEntry
 * providing access to all smoothed as well as raw values.
 */
public class EquantFeatureSmoother extends FeatureProvider {
    private final int windowSize;

    public EquantFeatureSmoother(int windowSize) {
        this.windowSize = windowSize;
    }

    @Override
    protected void processInternally(Structure structure) {

    }

    public static class EquantFeatureContainer extends FeatureContainerEntry {
        public EquantFeatureContainer(FeatureProvider featureProvider) {
            super(featureProvider);
        }
    }
}
