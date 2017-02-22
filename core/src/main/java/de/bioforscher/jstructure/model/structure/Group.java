package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.CIFParser;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * Created by S on 06.10.2016.
 */
public class Group extends AbstractFeatureContainer implements AtomContainer {
    /**
     * reference to an undefined group - this is used by atoms without explicit parent reference
     */
    static final Group UNKNOWN_GROUP = new Group("UNK", 0, GroupInformation.UNKNOWN_GROUP, false);

    private int residueNumber;
    private List<Atom> atoms;
    /**
     * Handle to the container element.
     */
    private Chain parentChain;
    private String identifier;
    private GroupInformation groupInformation;
    private boolean parentChainIsTermianted;

    public Group(String pdbName, int residueNumber, GroupInformation groupInformation, boolean parentChainIsTerminated) {
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        this.groupInformation = groupInformation;
        this.parentChainIsTermianted = parentChainIsTerminated;
    }

    Group() {

    }

    @Deprecated
    public Group(String pdbName, int residueNumber) {
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        //TODO it is not really nice, that the data model is actually 'parsing' stuff
        this.groupInformation = CIFParser.parseLigandInformation(pdbName);
    }

    public Group(Group group) {
        this.residueNumber = group.residueNumber;
        // deep clone entries
        this.atoms = group.atoms()
                .map(Atom::new)
                .collect(Collectors.toList());
        this.atoms().forEach(atom -> atom.setParentGroup(this));
        // reference parent
        this.parentChain = group.parentChain;
        this.groupInformation = group.groupInformation;
        // set reference to feature map
        setFeatureMap(group.getFeatureMap());
    }

    public GroupInformation getGroupInformation() {
        return groupInformation;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public Group(List<Atom> atoms) {
        this.atoms = atoms;
        this.groupInformation = GroupInformation.UNKNOWN_AMINO_ACID;
    }

    public boolean isAminoAcid() {
        return !parentChainIsTermianted && groupInformation.getType().contains("PEPTIDE LINKING");
    }

    public boolean isNucleotide() {
        return !parentChainIsTermianted && groupInformation.getType().contains("NA LINKING");
    }

    public boolean isLigand() {
        return !isAminoAcid() && !isNucleotide();
    }

    public List<Atom> getAtoms() {
        return atoms;
    }

    public String getThreeLetterCode() {
        return groupInformation.getThreeLetterCode();
    }

    public void setGroupInformation(GroupInformation groupInformation) {
        this.groupInformation = groupInformation;
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
        return parentChain != null ? parentChain : Chain.UNKNOWN_CHAIN;
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
        return identifier == null ? getThreeLetterCode() + "-" + residueNumber : identifier;
    }

    boolean isInTerminatedParentChain() {
        return parentChainIsTermianted;
    }
}