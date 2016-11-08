package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;

import java.util.*;

/**
 *
 * Created by S on 06.10.2016.
 */
public abstract class Group implements AtomContainer, AtomRecordWriter {
    protected List<Atom> atoms;
    protected Map<Enum, Object> featureMap;
    protected int residueNumber;
    /**
     * 3-letter pdb name
     */
    protected String pdbName;
    /**
     * Handle to the container element.
     */
    protected Chain parentChain;

    public Group(String pdbName, int residueNumber) {
        this.pdbName = pdbName;
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        this.featureMap = new HashMap<>();
    }

    @Override
    public List<Atom> getAtoms() {
        return atoms;
    }

    public String getPdbName() {
        return pdbName;
    }

    /**
     * Returns the {@link Chain}-wide unique number to identify this getResidue.
     * @return an integer
     */
    public int getResidueNumber() {
        return residueNumber;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " name='" + pdbName + "' resNum='" + residueNumber + "' size='" +
                atoms.size() + "'";
    }

    /**
     * Registers a child. This object will assign a reference to itself to the atom.
     * @param atom the atom to process
     */
    public void addAtom(Atom atom) {
        atoms.add(atom);
        atom.setParentGroup(this);
    }

    Optional<Atom> tryToGetAtomByName(final AtomNameFilter filter) {
        return atoms().filter(filter).findFirst();
    }

    Atom getAtomByName(final AtomNameFilter filter) {
        return tryToGetAtomByName(filter).orElseThrow(() -> new NoSuchElementException(String.format("no atom matching " +
                "%s in residue %s", filter.getAcceptedAtomNames(), toString())));
    }

    /**
     * Returns the {@link Chain} this getResidue is associated to.
     * @return the parent container
     */
    public Chain getParentChain() {
        return parentChain;
    }

    /**
     * Package-private method to set the parent reference.
     * @param parentChain the parent
     */
    void setParentChain(Chain parentChain) {
        this.parentChain = parentChain;
    }

    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}