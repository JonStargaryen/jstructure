package de.bioforscher.jstructure.model.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
    public void clear() {
        atoms.clear();
    }

    @Override
    public void clearPseudoAtoms() {
        atoms.removeIf(Atom::isVirtual);
    }

    public String getPdbName() {
        return pdbName;
    }

    /**
     * Returns the {@link Chain}-wide unique number to identify this residue.
     * @return an integer
     */
    public int getResidueNumber() {
        return residueNumber;
    }

    @Override
    public String composePDBRecord() {
        return atoms().map(Atom::composePDBRecord)
                      .collect(Collectors.joining(System.lineSeparator()));
    }

    /**
     * Registers a child. This object will assign a reference to itself to the atom.
     * @param atom the atom to process
     */
    public void addAtom(Atom atom) {
        atoms.add(atom);
        atom.setParentGroup(this);
    }

    /**
     * Returns the {@link Chain} this residue is associated to.
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

    @Override
    public Stream<Atom> atoms() {
        return atoms.stream();
    }

    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}