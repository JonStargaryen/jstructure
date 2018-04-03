package de.bioforscher.jstructure.efr.model;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

public class HotSpotScoring extends FeatureContainerEntry {
    private final int ecCount;
    private final double cumStrength;
    private final double ecStrength;
    private final int conservation;

    public HotSpotScoring(int ecCount, double cumStrength, double ecStrength, int conservation) {
        super(null);
        this.ecCount = ecCount;
        this.cumStrength = cumStrength;
        this.ecStrength = ecStrength;
        this.conservation = conservation;
    }

    public HotSpotScoring() {
        this(0, 0, 0, 0);
    }

    public int getEcCount() {
        return ecCount;
    }

    public double getCumStrength() {
        return cumStrength;
    }

    public double getEcStrength() {
        return ecStrength;
    }

    public int getConservation() {
        return conservation;
    }
}