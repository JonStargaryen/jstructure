package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

import java.util.List;

public class ResidueTopologicProperties extends FeatureContainerEntry {
    private final double betweenness;
    private final double closeness;
    private final double clusteringCoefficient;
    private final double distinctNeighborhoodCount;

    public ResidueTopologicProperties(FeatureProvider featureProvider,
                                      double betweenness,
                                      double closeness,
                                      double clusteringCoefficient,
                                      double distinctNeighborhoodCount) {
        super(featureProvider);
        this.betweenness = betweenness;
        this.closeness = closeness;
        this.clusteringCoefficient = clusteringCoefficient;
        this.distinctNeighborhoodCount = distinctNeighborhoodCount;
    }

    public ResidueTopologicProperties(List<ResidueTopologicProperties> residueTopologicProperties) {
        this(residueTopologicProperties.get(0).getFeatureProvider(),
                residueTopologicProperties.stream()
                        .mapToDouble(ResidueTopologicProperties::getBetweenness)
                        .average()
                        .getAsDouble(),
                residueTopologicProperties.stream()
                        .mapToDouble(ResidueTopologicProperties::getCloseness)
                        .average()
                        .getAsDouble(),
                residueTopologicProperties.stream()
                        .mapToDouble(ResidueTopologicProperties::getClusteringCoefficient)
                        .average()
                        .getAsDouble(),
                residueTopologicProperties.stream()
                        .mapToDouble(ResidueTopologicProperties::getDistinctNeighborhoodCount)
                        .average()
                        .getAsDouble());
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

    public double getDistinctNeighborhoodCount() {
        return distinctNeighborhoodCount;
    }
}
