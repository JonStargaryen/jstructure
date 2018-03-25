package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

import java.util.List;

@DefaultFeatureProvider(TopologicPropertyCalculator.class)
public class ResidueTopologicPropertiesContainer extends FeatureContainerEntry {
    private final ResidueTopologicProperties fullPlip;
    private final ResidueTopologicProperties hydrogenPlip;
    private final ResidueTopologicProperties hydrophobicPlip;
    private final ResidueTopologicProperties conventional;

    ResidueTopologicPropertiesContainer(FeatureProvider featureProvider,
                                        ResidueTopologicProperties fullPlip,
                                        ResidueTopologicProperties hydrogenPlip,
                                        ResidueTopologicProperties hydrophobicPlip,
                                        ResidueTopologicProperties conventional) {
        super(featureProvider);
        this.fullPlip = fullPlip;
        this.hydrogenPlip = hydrogenPlip;
        this.hydrophobicPlip = hydrophobicPlip;
        this.conventional = conventional;
    }

    public ResidueTopologicProperties getFullPlip() {
        return fullPlip;
    }

    public ResidueTopologicProperties getHydrogenPlip() {
        return hydrogenPlip;
    }

    public ResidueTopologicProperties getHydrophobicPlip() {
        return hydrophobicPlip;
    }

    public ResidueTopologicProperties getConventional() {
        return conventional;
    }

    public static class ResidueTopologicProperties extends FeatureContainerEntry {
        private final double betweenness;
        private final double closeness;
        private final double clusteringCoefficient;
        private final double distinctNeighborhoodCount;

        ResidueTopologicProperties(FeatureProvider featureProvider,
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
                            .orElse(0.0),
                    residueTopologicProperties.stream()
                            .mapToDouble(ResidueTopologicProperties::getCloseness)
                            .average()
                            .orElse(0.0),
                    residueTopologicProperties.stream()
                            .mapToDouble(ResidueTopologicProperties::getClusteringCoefficient)
                            .average()
                            .orElse(0.0),
                    residueTopologicProperties.stream()
                            .mapToDouble(ResidueTopologicProperties::getDistinctNeighborhoodCount)
                            .average()
                            .orElse(0.0));
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
}
