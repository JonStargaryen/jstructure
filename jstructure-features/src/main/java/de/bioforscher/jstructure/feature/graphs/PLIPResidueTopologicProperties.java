package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureProvider;

public class PLIPResidueTopologicProperties extends ResidueTopologicProperties {
    public PLIPResidueTopologicProperties(FeatureProvider featureProvider,
                                          double betweenness,
                                          double closeness,
                                          double clusteringCoefficient) {
        super(featureProvider, betweenness, closeness, clusteringCoefficient);
    }
}
