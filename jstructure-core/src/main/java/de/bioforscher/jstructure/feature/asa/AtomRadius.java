package de.bioforscher.jstructure.feature.asa;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The radius of an atom, used to compute the accessible surface area.
 * Created by bittrich on 5/17/17.
 */
public class AtomRadius extends FeatureContainerEntry {
    private final double radius;

    public AtomRadius(AbstractFeatureProvider featureProvider, double radius) {
        super(featureProvider);
        this.radius = radius;
    }

    public double getRadius() {
        return radius;
    }
}
