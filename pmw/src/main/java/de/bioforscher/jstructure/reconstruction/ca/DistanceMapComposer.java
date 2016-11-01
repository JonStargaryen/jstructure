package de.bioforscher.jstructure.reconstruction.ca;

import de.bioforscher.jstructure.model.structure.Protein;

/**
 * Converts an existing protein into its distance map representation.
 * Created by S on 26.10.2016.
 */
public class DistanceMapComposer {
    public enum FeatureNames {
        DISTANCE_MAP
    }

    public double[][] computeDistanceMap(Protein protein) {
        double[][] distanceMap = new double[protein.getSize()][protein.getSize()];
        protein.residuePairs();
//        for(Chain chain1 : protein.chains) {
//            for(Residue residue1 : chain1.residues) {
//                for(Chain chain2 : protein.chains) {
//                    for(Residue residue2 : chain2.residues) {
//                        double distance = this.linearAlgebra.distance(this.modelConverter.getCA(residue1).xyz,
//                                this.modelConverter.getCA(residue2).xyz);
//                        int i = residue1.residueId;
//                        int j = residue2.residueId;
//                        distanceMap[i][j] = distance;
//                        distanceMap[j][i] = distance;
//                    }
//                }
//            }
//        }
        return distanceMap;
    }
}