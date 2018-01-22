package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

public class ResidueTopologicProperties extends FeatureContainerEntry {
    private final double betweenness;
    private final double closeness;
    private final double clusteringCoefficient;

    public ResidueTopologicProperties(FeatureProvider featureProvider,
                                      double betweenness,
                                      double closeness,
                                      double clusteringCoefficient) {
        super(featureProvider);
        this.betweenness = betweenness;
        this.closeness = closeness;
        this.clusteringCoefficient = clusteringCoefficient;
    }

    public double getBetweenness() {
        return betweenness;
    }

    public double getCloseness() {
        return closeness;
    }

    public double getClusteringCoefficient() {
        return clusteringCoefficient;
    }
}
