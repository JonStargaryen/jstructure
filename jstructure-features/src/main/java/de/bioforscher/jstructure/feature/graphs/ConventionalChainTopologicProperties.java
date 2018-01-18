package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureProvider;

public class ConventionalChainTopologicProperties extends ChainTopologicProperties {
    public ConventionalChainTopologicProperties(FeatureProvider featureProvider, double averagePathLength, double clusteringCoefficient) {
        super(featureProvider, averagePathLength, clusteringCoefficient);
    }
}
