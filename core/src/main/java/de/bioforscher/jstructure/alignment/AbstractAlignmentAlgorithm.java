package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;
import java.util.Objects;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    private List<String> alignedAtomNames;

    protected AbstractAlignmentAlgorithm(List<String> alignedAtomNames) {
        this.alignedAtomNames = alignedAtomNames;
    }

    public AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        return alignInternal(prepareAtomContainer(atomContainer1), prepareAtomContainer(atomContainer2));
    }

    protected abstract AlignmentResult alignInternal(AtomContainer atomContainer1, AtomContainer atomContainer2);

    /**
     * Selects a subset of this container containing only atoms which will actually be align (specified by
     * {@link #getAlignedAtomNames()}). Furthermore, ensures correct ordering of the atoms.
     * @return the subset of atoms to align
     */
    private AtomContainer prepareAtomContainer(AtomContainer atomContainer) {
        // sort atoms in each residue container
        atomContainer.atoms()
                // move to parent group level to be able to sort atoms by name
                .map(Atom::getParentGroup)
                .filter(Objects::nonNull)
                // multiple atoms will point to the same parent
                .distinct()
                .filter(Group::isAminoAcid)
                .forEach(AminoAcid::orderAtomsByName);
        // select, clone and wrap in container
        return Selection.on(atomContainer)
                        .atomName(alignedAtomNames)
                        // request the filtered atoms to be cloned, so we can manipulate there coordinates safely
                        .cloneElements()
                        .asAtomContainer();
    }

    @Override
    public List<String> getAlignedAtomNames() {
        return alignedAtomNames;
    }
}
