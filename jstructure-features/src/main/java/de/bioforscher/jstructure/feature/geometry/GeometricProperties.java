package de.bioforscher.jstructure.feature.geometry;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

@DefaultFeatureProvider(GeometricPropertyCalculator.class)
public class GeometricProperties extends FeatureContainerEntry {
    private final double distanceToCenterOfMass;
    private final double distanceToCentroid;

    public GeometricProperties(FeatureProvider featureProvider,
                               double distanceToCenterOfMass,
                               double distanceToCentroid) {
        super(featureProvider);
        this.distanceToCenterOfMass = distanceToCenterOfMass;
        this.distanceToCentroid = distanceToCentroid;
    }

    public double getDistanceToCenterOfMass() {
        return distanceToCenterOfMass;
    }

    public double getDistanceToCentroid() {
        return distanceToCentroid;
    }
}
