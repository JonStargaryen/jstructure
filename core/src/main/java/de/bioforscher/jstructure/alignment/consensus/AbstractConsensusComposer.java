package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;

/**
 * The abstract capabilities of algorithms which merge several observations of fragments or motifs into 1 consensus
 * representation.
 * Created by S on 05.12.2016.
 */
@Deprecated
public abstract class AbstractConsensusComposer {
    /**
     * Merges two atoms by updating the coordinates of the first argument.
     * @param atomToManipulate this atom's coordinates will be updated
     * @param atomToConsume this atom's coordinates will be merged upon the first argument
     */
    private static void mergeAtomsByCentroid(Atom atomToManipulate, Atom atomToConsume) {
        atomToManipulate.setCoordinates(LinearAlgebra3D.divide(LinearAlgebra3D.add(atomToManipulate.getCoordinates(), atomToConsume.getCoordinates()), 2));
    }

    /**
     * Merges two groups (i.e. merging their respective atoms) by updating the coordinates of the first argument.
     * @param groupToManipulate this group's coordinates will be updated
     * @param groupToConsume this group's coordinates will be merged upon the first argument
     */
    private static void mergeGroupsByCentroid(Group groupToManipulate, Group groupToConsume) {
        for(int atomIndex = 0; atomIndex < groupToManipulate.getAtoms().size(); atomIndex++) {
            mergeAtomsByCentroid(groupToManipulate.getAtoms().get(atomIndex), groupToConsume.getAtoms().get(atomIndex));
        }
    }

    /**
     * Merges two containers by determining the shared set of atoms and updating their coordinates.
     * @param reference the container whose atom and group names and indices will be kept - the container will not be
     *                  manipulated in any way however
     * @param containerToConsume the container to merge upon the first
     * @return a distinct clone of the first argument with updated coordinates (containing only the subset of compatible
     * atoms)
     */
    public static GroupContainer mergeContainersByCentroid(AtomContainer reference, AtomContainer containerToConsume) {
        // get intersecting pairs of atoms wrapped in an atom container
        // 12/14/16 - moved to backboneOnly flag
        Pair<GroupContainer, GroupContainer> groupContainerPair =
                LinearAlgebraAtom.comparableGroupContainerPair(reference,
                        containerToConsume,
                        AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES,
                        AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES);

        for(int groupIndex = 0; groupIndex < groupContainerPair.getLeft().getGroups().size(); groupIndex++) {
            mergeGroupsByCentroid(groupContainerPair.getLeft().getGroups().get(groupIndex),
                    groupContainerPair.getRight().getGroups().get(groupIndex));
        }
        return groupContainerPair.getLeft();
    }
}
