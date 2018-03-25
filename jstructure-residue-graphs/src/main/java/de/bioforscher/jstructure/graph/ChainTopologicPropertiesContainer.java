package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

@DefaultFeatureProvider(TopologicPropertyCalculator.class)
public class ChainTopologicPropertiesContainer extends FeatureContainerEntry {
    private final ChainTopologicProperties fullPlip;
    private final ChainTopologicProperties hydrogenPlip;
    private final ChainTopologicProperties hydrophobicPlip;
    private final ChainTopologicProperties conventional;

    ChainTopologicPropertiesContainer(FeatureProvider featureProvider,
                                      ChainTopologicProperties fullPlip,
                                      ChainTopologicProperties hydrogenPlip,
                                      ChainTopologicProperties hydrophobicPlip,
                                      ChainTopologicProperties conventional) {
        super(featureProvider);
        this.fullPlip = fullPlip;
        this.hydrogenPlip = hydrogenPlip;
        this.hydrophobicPlip = hydrophobicPlip;
        this.conventional = conventional;
    }

    public ChainTopologicProperties getFullPlip() {
        return fullPlip;
    }

    public ChainTopologicProperties getHydrogenPlip() {
        return hydrogenPlip;
    }

    public ChainTopologicProperties getHydrophobicPlip() {
        return hydrophobicPlip;
    }

    public ChainTopologicProperties getConventional() {
        return conventional;
    }

    public static class ChainTopologicProperties extends FeatureContainerEntry {
        private final double averagePathLength;
        private final double clusteringCoefficient;

        ChainTopologicProperties(FeatureProvider featureProvider,
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

}
