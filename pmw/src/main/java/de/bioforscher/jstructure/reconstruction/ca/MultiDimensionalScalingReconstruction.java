package de.bioforscher.jstructure.reconstruction.ca;

import de.bioforscher.jstructure.mathematics.mds.MultiDimensionalScaling;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
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
        List<Residue> residues = protein.residues().collect(Collectors.toList());

        // create mapping using MDS
        MultiDimensionalScaling mds = new MultiDimensionalScaling();
        double[][] distanceMap = protein.getFeature(double[][].class, DistanceMapComposer.FeatureNames.DISTANCE_MAP);
        List<double[]> placedAtoms = mds.computeEmbedding(distanceMap);

        Pair.sequentialPairsOf(residues, placedAtoms)
            .forEach(pair -> {
                Residue residue = pair.getLeft();
                residue.addAtom(new Atom("CA", 0, Element.C, pair.getRight()));
            });
    }
}
