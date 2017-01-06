package de.bioforscher.jstructure.reconstruction.ca;

import de.bioforscher.jstructure.mathematics.mds.MultiDimensionalScaling;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.reconstruction.ReconstructionAlgorithm;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A rather basic reconstruction algorithm for the placement of CA atoms by knowing their atomic distances. Finds a
 * configuration of 3D coordinates explaining the given distances as well as possible.
 * Created by S on 26.10.2016.
 */
public class MultiDimensionalScalingReconstruction implements ReconstructionAlgorithm {
    @Override
    public void reconstruct(Protein protein) {
        // ensure container is empty
        protein.clearAtoms();
        List<Group> residues = Selection.on(protein)
                .aminoAcids()
                .cloneElements()
                .asFilteredGroups()
                .collect(Collectors.toList());

        // create mapping using MDS
        MultiDimensionalScaling mds = new MultiDimensionalScaling();
        double[][] distanceMap = protein.getFeature(double[][].class, DistanceMapComposer.FeatureNames.DISTANCE_MAP);
        List<double[]> placedAtoms = mds.computeEmbedding(distanceMap);

        IntStream.range(0, residues.size())
                .forEach(i -> residues.get(i).addAtom(new Atom(AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME, 0, Element.C, placedAtoms.get(i))));
    }
}
