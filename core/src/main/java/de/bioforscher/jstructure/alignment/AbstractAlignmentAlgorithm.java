package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;

import java.util.Objects;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    private AtomNameFilter alignedAtomNameFilter;

    AbstractAlignmentAlgorithm(AtomNameFilter alignedAtomNameFilter) {
        this.alignedAtomNameFilter = alignedAtomNameFilter;
    }

    public AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        return alignInternal(prepareAtomContainer(atomContainer1), prepareAtomContainer(atomContainer2));
    }

    abstract AlignmentResult alignInternal(AtomContainer atomContainer1, AtomContainer atomContainer2);

    /**
     * Selects a subset of this container containing only atoms which will actually be align (specified by
     * {@link #getAlignedAtomNameFilter()}). Furthermore, ensures correct ordering of the atoms.
     * @return the subset of atoms to align
     */
    private AtomContainer prepareAtomContainer(AtomContainer atomContainer) {
        // sort atoms in each residue container
        atomContainer.atoms()
                .map(Atom::getParentGroup)
                .filter(Objects::nonNull)
                .distinct()
                .map(Residue.class::cast)
                .forEach(AminoAcid::orderAtomsByName);
        // filter and wrap in container
        return AtomContainer.of(atomContainer.atoms()
                //TODO clone here?
                .filter(getAlignedAtomNameFilter()));
    }

    @Override
    public AtomNameFilter getAlignedAtomNameFilter() {
        return alignedAtomNameFilter;
    }
}
