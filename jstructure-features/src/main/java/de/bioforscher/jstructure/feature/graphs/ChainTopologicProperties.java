package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

public abstract class ChainTopologicProperties extends FeatureContainerEntry {
    private final double averagePathLength;
    private final double clusteringCoefficient;

    public ChainTopologicProperties(FeatureProvider featureProvider,
                                    double averagePathLength,
                                    double clusteringCoefficient) {
        super(featureProvider);
        this.averagePathLength = averagePathLength;
        this.clusteringCoefficient = clusteringCoefficient;
    }

    public double getAveragePathLength() {
        return averagePathLength;
    }

    public double getClusteringCoefficient() {
        return clusteringCoefficient;
    }
}
