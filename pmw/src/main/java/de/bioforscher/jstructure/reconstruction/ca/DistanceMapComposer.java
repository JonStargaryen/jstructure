package de.bioforscher.jstructure.reconstruction.ca;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.List;
import java.util.stream.Collectors;

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
        // we need a mapping from residues to their index in the containing list
        List<Residue> residues = protein.residues().collect(Collectors.toList());

        protein.residuePairs()
               .forEach(pair -> {
                   Residue first = pair.getLeft();
                   Residue second = pair.getRight();
                   int indexFirst = residues.indexOf(first);
                   int indexSecond = residues.indexOf(second);
                   distanceMap[indexFirst][indexSecond] = LinearAlgebra3D.distance(first.getAlphaCarbon().getCoordinates(),
                           second.getAlphaCarbon().getCoordinates());
               });

        return distanceMap;
    }
}