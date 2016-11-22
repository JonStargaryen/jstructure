package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * Created by S on 06.10.2016.
 */
public class Group extends AbstractFeatureContainer implements AtomContainer {
    public enum GroupType {
        AMINO_ACID,
        NUCLEOTIDE,
        HETATM
    }
    private int residueNumber;
    private List<Atom> atoms;
    /**
     * 3-letter pdb name
     */
    private String pdbName;
    /**
     * Handle to the container element.
     */
    private Chain parentChain;
    private GroupType groupType;
    private String identifier;

    public Group(String pdbName, int residueNumber) {
        this.pdbName = pdbName;
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        this.groupType = determineGroupType();
    }

    public Group(Group group) {
        this.residueNumber = group.residueNumber;
        // deep clone entries
        this.atoms = group.atoms()
                .map(Atom::new)
                .collect(Collectors.toList());
        this.pdbName = group.pdbName;
        // reference parent
        this.parentChain = group.parentChain;
        this.groupType = group.groupType;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public Group(List<Atom> atoms) {
        this.atoms = atoms;
    }

    private GroupType determineGroupType() {
        //TODO impl
        return GroupType.AMINO_ACID;
    }

    public boolean isAminoAcid() {
        return groupType.equals(GroupType.AMINO_ACID);
    }

    public boolean isNucleotide() {
        return groupType.equals(GroupType.NUCLEOTIDE);
    }

    public boolean isHetatm() {
        return groupType.equals(GroupType.HETATM);
    }

    public List<Atom> getAtoms() {
        return atoms;
    }

    public String getPdbName() {
        return pdbName;
    }

    public void setPdbName(String pdbName) {
        this.pdbName = pdbName;
    }

    public GroupType getGroupType() {
        return groupType;
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
        return getClass().getSimpleName() + " identifier='" + getIdentifier() + "' size='" + getAtoms().size() + "'";
    }

    /**
     * Registers a child. This object will assign a reference to itself to the atom.
     * @param atom the atom to process
     */
    public void addAtom(Atom atom) {
        getAtoms().add(atom);
        atom.setParentGroup(this);
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

    @Override
    public String getIdentifier() {
        return identifier == null ? pdbName + "-" + residueNumber : identifier;
    }
}