package de.bioforscher.jstructure.reconstruction.ca;

import de.bioforscher.jstructure.mathematics.mds.MultiDimensionalScaling;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.reconstruction.ReconstructionAlgorithm;

import java.util.List;
import java.util.stream.Collectors;

/**
 * A rather basic reconstruction algorithm for the placement of CA atoms by knowing their atomic distances. Finds a
 * configuration of 3D coordinates explaining the given distances as well as possible.
 * Created by S on 26.10.2016.
 */
public class MultiDimensionalScalingReconstruction implements ReconstructionAlgorithm {
    @Override
    public void reconstruct(Protein protein) {
        // ensure container is empty
        protein.clear();

        MultiDimensionalScaling mds = new MultiDimensionalScaling();
        double[][] distanceMap = protein.getFeature(double[][].class, DistanceMapComposer.FeatureNames.DISTANCE_MAP);
        List<double[]> placedAtoms = mds.computeEmbedding(distanceMap);

        Pair.sequentialPairsOf(protein.residues().collect(Collectors.toList()), placedAtoms);
//        for(int i = 0; i < placedAtoms.size(); i++) {
//            this.modelConverter.createAtom(residues.get(i), ModelConverter.BACKBONE_CA_NAME, placedAtoms.get(i));
//        }
    }
}
