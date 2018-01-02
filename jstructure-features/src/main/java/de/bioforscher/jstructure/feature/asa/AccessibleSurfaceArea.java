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
    /**
     * see Rost B, Sander C. Conservation and prediction of solvent accessibility in protein families. Proteins 1994 Nov;20(3):216-26.
     */
    private static final double BURIED_RASA_THRESHOLD = 0.16;
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

    public boolean isBuried() {
        return relativeAccessibleSurfaceArea < BURIED_RASA_THRESHOLD;
    }

    public boolean isExposed() {
        return !isBuried();
    }
}
