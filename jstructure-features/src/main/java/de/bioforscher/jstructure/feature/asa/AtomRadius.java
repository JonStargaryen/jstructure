package de.bioforscher.jstructure.feature.asa;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

/**
 * The radius of an atom, used to calculate the accessible surface area.
 * Created by bittrich on 5/17/17.
 */
@DefaultFeatureProvider(AccessibleSurfaceAreaCalculator.class)
public class AtomRadius extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double radius;

    AtomRadius(FeatureProvider featureProvider, double radius) {
        super(featureProvider);
        this.radius = radius;
    }

    public double getRadius() {
        return radius;
    }

    @Override
    public Double getValue() {
        return getRadius();
    }
}
