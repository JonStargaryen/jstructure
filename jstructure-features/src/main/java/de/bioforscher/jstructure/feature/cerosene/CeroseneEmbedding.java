package de.bioforscher.jstructure.feature.cerosene;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

/**
 * The sequence cerosene embedding of groups.
 * Created by bittrich on 5/17/17.
 */
@DefaultFeatureProvider(SequenceCerosene.class)
public class CeroseneEmbedding extends FeatureContainerEntry {
    private final double[] rgb;
    private final double[] hsv;

    CeroseneEmbedding(FeatureProvider featureProvider, double[] rgb, double[] hsv) {
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
