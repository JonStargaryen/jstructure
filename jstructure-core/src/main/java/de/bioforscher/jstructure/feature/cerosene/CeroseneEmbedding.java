package de.bioforscher.jstructure.feature.cerosene;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

/**
 * The sequence cerosene embedding of groups.
 * Created by bittrich on 5/17/17.
 */
public class CeroseneEmbedding extends FeatureContainerEntry {
    private final double[] rgb;
    private final double[] hsv;

    public CeroseneEmbedding(AbstractFeatureProvider featureProvider, double[] rgb, double[] hsv) {
        super(featureProvider);
        this.rgb = rgb;
        this.hsv = hsv;
    }

    public double[] getRgb() {
        return rgb;
    }

    public double[] getHsv() {
        return hsv;
    }
}
