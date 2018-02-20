package de.bioforscher.start2fold.model.vector;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.ArrayList;
import java.util.List;

public class SmoothedFeatureVector extends FeatureVector {
    public SmoothedFeatureVector(double secondaryStructureElementSize, double localHydrogen, double localHydrophobic, double localBackbone, double localInteractions, double nonLocalHydrogen, double nonLocalHydrophobic, double nonLocalBackbone, double nonLocalInteractions, double energy, double egor, double rasa, double betweenness, double closeness, double clusteringCoefficient, double hydrogenBetweenness, double hydrogenCloseness, double hydrogenClusteringCoefficient, double hydrophobicBetweenness, double hydrophobicCloseness, double hydrophobicClusteringCoefficient, double convBetweenness, double convCloseness, double convClusteringCoefficient, double distinctNeighborhoods, double convDistinctNeighborhoods) {
        super(secondaryStructureElementSize, localHydrogen, localHydrophobic, localBackbone, localInteractions, nonLocalHydrogen, nonLocalHydrophobic, nonLocalBackbone, nonLocalInteractions, energy, egor, rasa, betweenness, closeness, clusteringCoefficient, hydrogenBetweenness, hydrogenCloseness, hydrogenClusteringCoefficient, hydrophobicBetweenness, hydrophobicCloseness, hydrophobicClusteringCoefficient, convBetweenness, convCloseness, convClusteringCoefficient, distinctNeighborhoods, convDistinctNeighborhoods);
    }

    public static void assignSmoothedFeatureVector(List<AminoAcid> aminoAcids, AminoAcid aminoAcid) {
        int index = aminoAcid.getResidueIndex();
        List<AminoAcid> aminoAcidsToSmooth = new ArrayList<>();
        aminoAcidsToSmooth.add(aminoAcid);

        // get N-terminal residues
        for(int i = 0; i < 4; i++) {
            int indexToGet = index - (i + 1);
            if(indexToGet > -1) {
                aminoAcidsToSmooth.add(aminoAcids.get(indexToGet));
            }
        }

        // get C-terminal residues
        for(int i = 0; i < 4; i++) {
            int indexToGet = index + (i + 1);
            if(indexToGet < aminoAcids.size()) {
                aminoAcidsToSmooth.add(aminoAcids.get(indexToGet));
            }
        }

        // assign smoothed values
        ValueSmoother.smoothValues(aminoAcidsToSmooth, aminoAcid);
    }
}