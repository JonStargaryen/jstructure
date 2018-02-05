package de.bioforscher.start2fold.model.vector;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

public abstract class FeatureVector extends FeatureContainerEntry {
    private final double secondaryStructureElementSize;
    private final double localHydrogen;
    private final double localHydrophobic;
    private final double localBackbone;
    private final double localInteractions;
    private final double nonLocalHydrogen;
    private final double nonLocalHydrophobic;
    private final double nonLocalBackbone;
    private final double nonLocalInteractions;
    private final double energy;
    private final double egor;
    //        private final double eccount;
//        private final double cumstrength;
//        private final double ecstrength;
//        private final double conservation;
    private final double rasa;
    private final double betweenness;
    private final double closeness;
    private final double clusteringCoefficient;
    private final double hydrogenBetweenness;
    private final double hydrogenCloseness;
    private final double hydrogenClusteringCoefficient;
    private final double hydrophobicBetweenness;
    private final double hydrophobicCloseness;
    private final double hydrophobicClusteringCoefficient;
    private final double convBetweenness;
    private final double convCloseness;
    private final double convClusteringCoefficient;
    private final double distinctNeighborhoods;
    private final double convDistinctNeighborhoods;

    FeatureVector(double secondaryStructureElementSize,
                  double localHydrogen,
                  double localHydrophobic,
                  double localBackbone,
                  double localInteractions,
                  double nonLocalHydrogen,
                  double nonLocalHydrophobic,
                  double nonLocalBackbone,
                  double nonLocalInteractions,
                  double energy,
                  double egor,
//                             double eccount,
//                             double cumstrength,
//                             double ecstrength,
//                             double conservation,
                  double rasa,
                  double betweenness,
                  double closeness,
                  double clusteringCoefficient,
                  double hydrogenBetweenness,
                  double hydrogenCloseness,
                  double hydrogenClusteringCoefficient,
                  double hydrophobicBetweenness,
                  double hydrophobicCloseness,
                  double hydrophobicClusteringCoefficient,
                  double convBetweenness,
                  double convCloseness,
                  double convClusteringCoefficient,
                  double distinctNeighborhoods,
                  double convDistinctNeighborhoods) {
        super(null);
        this.secondaryStructureElementSize = secondaryStructureElementSize;
        this.localHydrogen = localHydrogen;
        this.localHydrophobic = localHydrophobic;
        this.localBackbone = localBackbone;
        this.localInteractions = localInteractions;
        this.nonLocalHydrogen = nonLocalHydrogen;
        this.nonLocalHydrophobic = nonLocalHydrophobic;
        this.nonLocalBackbone = nonLocalBackbone;
        this.nonLocalInteractions = nonLocalInteractions;
        this.energy = energy;
        this.egor = egor;
//            this.eccount = eccount;
//            this.cumstrength = cumstrength;
//            this.ecstrength = ecstrength;
//            this.conservation = conservation;
        this.rasa = rasa;
        this.betweenness = betweenness;
        this.closeness = closeness;
        this.clusteringCoefficient = clusteringCoefficient;
        this.hydrogenBetweenness = hydrogenBetweenness;
        this.hydrogenCloseness = hydrogenCloseness;
        this.hydrogenClusteringCoefficient = hydrogenClusteringCoefficient;
        this.hydrophobicBetweenness = hydrophobicBetweenness;
        this.hydrophobicCloseness = hydrophobicCloseness;
        this.hydrophobicClusteringCoefficient = hydrophobicClusteringCoefficient;
        this.convBetweenness = convBetweenness;
        this.convCloseness = convCloseness;
        this.convClusteringCoefficient = convClusteringCoefficient;
        this.distinctNeighborhoods = distinctNeighborhoods;
        this.convDistinctNeighborhoods = convDistinctNeighborhoods;
    }

    public double getConvBetweenness() {
        return convBetweenness;
    }

    public double getConvCloseness() {
        return convCloseness;
    }

    public double getConvClusteringCoefficient() {
        return convClusteringCoefficient;
    }

    public double getConvDistinctNeighborhoods() {
        return convDistinctNeighborhoods;
    }

    public double getSecondaryStructureElementSize() {
        return secondaryStructureElementSize;
    }

    public double getLocalHydrogen() {
        return localHydrogen;
    }

    public double getLocalHydrophobic() {
        return localHydrophobic;
    }

    public double getLocalBackbone() {
        return localBackbone;
    }

    public double getLocalInteractions() {
        return localInteractions;
    }

    public double getNonLocalHydrogen() {
        return nonLocalHydrogen;
    }

    public double getNonLocalHydrophobic() {
        return nonLocalHydrophobic;
    }

    public double getNonLocalBackbone() {
        return nonLocalBackbone;
    }

    public double getNonLocalInteractions() {
        return nonLocalInteractions;
    }

    public double getEnergy() {
        return energy;
    }

    public double getEgor() {
        return egor;
    }

//        public double getEccount() {
//            return eccount;
//        }
//
//        public double getCumstrength() {
//            return cumstrength;
//        }
//
//        public double getEcstrength() {
//            return ecstrength;
//        }
//
//        public double getConservation() {
//            return conservation;
//        }

    public double getRasa() {
        return rasa;
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

    public double getHydrogenBetweenness() {
        return hydrogenBetweenness;
    }

    public double getHydrogenCloseness() {
        return hydrogenCloseness;
    }

    public double getHydrogenClusteringCoefficient() {
        return hydrogenClusteringCoefficient;
    }

    public double getHydrophobicBetweenness() {
        return hydrophobicBetweenness;
    }

    public double getHydrophobicCloseness() {
        return hydrophobicCloseness;
    }

    public double getHydrophobicClusteringCoefficient() {
        return hydrophobicClusteringCoefficient;
    }

    public double getDistinctNeighborhoods() {
        return distinctNeighborhoods;
    }
}