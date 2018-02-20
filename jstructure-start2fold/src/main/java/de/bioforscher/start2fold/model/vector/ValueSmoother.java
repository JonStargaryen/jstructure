package de.bioforscher.start2fold.model.vector;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class ValueSmoother {
    public static void smoothValues(List<AminoAcid> aminoAcidsToSmooth, AminoAcid aminoAcid) {
        double secondaryStructureElementSize = smoothValue(aminoAcidsToSmooth, FeatureVector::getSecondaryStructureElementSize);

        double localHydrogen = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalHydrogen);
        double localHydrophobic = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalHydrophobic);
        double localBackbone = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalBackbone);
        double localInteractions = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalInteractions);

        double nonLocalHydrogen = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalHydrogen);
        double nonLocalHydrophobic = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalHydrophobic);
        double nonLocalBackbone = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalBackbone);
        double nonLocalInteractions = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalInteractions);

        double energy = smoothValue(aminoAcidsToSmooth, FeatureVector::getEnergy);
        double egor = smoothValue(aminoAcidsToSmooth, FeatureVector::getEgor);

        double rasa = smoothValue(aminoAcidsToSmooth, FeatureVector::getRasa);

        double betweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getBetweenness);
        double closeness = smoothValue(aminoAcidsToSmooth, FeatureVector::getCloseness);
        double clusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getClusteringCoefficient);

        double hydrogenBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenBetweenness);
        double hydrogenCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenCloseness);
        double hydrogenClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenClusteringCoefficient);

        double hydrophobicBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicBetweenness);
        double hydrophobicCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicCloseness);
        double hydrophobicClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicClusteringCoefficient);

        double convBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvBetweenness);
        double convCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvCloseness);
        double convClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvClusteringCoefficient);

        double distinctNeighborhoods = smoothValue(aminoAcidsToSmooth, FeatureVector::getDistinctNeighborhoods);
        double convDistinctNeighborhoods = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvDistinctNeighborhoods);

        SmoothedFeatureVector smoothedFeatureVector = new SmoothedFeatureVector(secondaryStructureElementSize,

                localHydrogen,
                localHydrophobic,
                localBackbone,
                localInteractions,

                nonLocalHydrogen,
                nonLocalHydrophobic,
                nonLocalBackbone,
                nonLocalInteractions,

                energy,
                egor,

                rasa,

                betweenness,
                closeness,
                clusteringCoefficient,

                hydrogenBetweenness,
                hydrogenCloseness,
                hydrogenClusteringCoefficient,

                hydrophobicBetweenness,
                hydrophobicCloseness,
                hydrophobicClusteringCoefficient,

                convBetweenness,
                convCloseness,
                convClusteringCoefficient,

                distinctNeighborhoods,
                convDistinctNeighborhoods);

        aminoAcid.getFeatureContainer().addFeature(smoothedFeatureVector);
    }

    private static double smoothValue(List<AminoAcid> aminoAcidsToSmooth, ToDoubleFunction<FeatureVector> mapping) {
        return aminoAcidsToSmooth.stream()
                .map(aminoAcid -> aminoAcid.getFeature(RawFeatureVector.class))
                .mapToDouble(mapping)
                .average()
                .orElseThrow(() -> new IllegalArgumentException("could not compute average"));
    }
}
