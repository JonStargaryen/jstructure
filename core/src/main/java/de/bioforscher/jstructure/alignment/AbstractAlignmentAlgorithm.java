package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    public AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        Pair<AtomContainer, AtomContainer> preparedPair = prepareAtomContainer(atomContainer1, atomContainer2);
        return alignInternal(preparedPair.getLeft(), preparedPair.getRight());
    }

    protected abstract AlignmentResult alignInternal(AtomContainer atomContainer1, AtomContainer atomContainer2);

    public static Pair<AtomContainer, AtomContainer> prepareAtomContainer(AtomContainer atomContainer1,
                                                                          AtomContainer atomContainer2) {

        return null;
    }

    private static GroupContainer cloneIntoGroupContainer(AtomContainer atomContainer) {
        return null;
    }

    private static void sortAtoms(AtomContainer atomContainer) {
        // sort atoms in each residue container
        atomContainer.atoms()
                // move to parent group level to be able to sort atoms by name
                .map(Atom::getParentGroup)
                // multiple atoms will point to the same parent
                .distinct()
                .forEach(AminoAcid::orderAtomsByName);
    }
}
