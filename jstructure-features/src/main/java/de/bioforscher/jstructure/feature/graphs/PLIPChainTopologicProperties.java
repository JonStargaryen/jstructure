package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureProvider;

public class PLIPChainTopologicProperties extends ChainTopologicProperties {
    public PLIPChainTopologicProperties(FeatureProvider featureProvider, double averagePathLength, double clusteringCoefficient) {
        super(featureProvider, averagePathLength, clusteringCoefficient);
    }
}
