package de.bioforscher.jstructure.feature.asa;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The entry for the accessible surface area.
 * Created by bittrich on 5/17/17.
 */
@DefaultFeatureProvider(AccessibleSurfaceAreaCalculator.class)
public class AccessibleSurfaceArea extends FeatureContainerEntry {
    private final double accessibleSurfaceArea;
    private final double relativeAccessibleSurfaceArea;

    AccessibleSurfaceArea(FeatureProvider featureProvider, double accessibleSurfaceArea, double relativeAccessibleSurfaceArea) {
        super(featureProvider);
        this.accessibleSurfaceArea = accessibleSurfaceArea;
        this.relativeAccessibleSurfaceArea = relativeAccessibleSurfaceArea;
    }

    public double getAccessibleSurfaceArea() {
        return accessibleSurfaceArea;
    }

    public double getRelativeAccessibleSurfaceArea() {
        return relativeAccessibleSurfaceArea;
    }
}
