package de.bioforscher.jstructure.membrane;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

public class Kink extends FeatureContainerEntry {
    private final double angle;
    private final boolean significantKink;

    public Kink(double angle, boolean significantKink) {
        super(null);
        this.angle = angle;
        this.significantKink = significantKink;
    }

    public double getAngle() {
        return angle;
    }

    public boolean isSignificantKink() {
        return significantKink;
    }
}