package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.CIFParser;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * Created by S on 06.10.2016.
 */
public class Group extends AbstractFeatureable implements AtomContainer {
    /**
     * reference to an undefined group - this is used by atoms without explicit parent reference
     */
    static final Group UNKNOWN_GROUP = new Group(0, GroupInformation.UNKNOWN_GROUP, false, false);

    private int residueNumber;
    private List<Atom> atoms;
    /**
     * Handle to the container element.
     */
    private Chain parentChain;
    private String identifier;
    private GroupInformation groupInformation;
    private boolean parentChainIsTerminated;
    private boolean hetAtm;

    public Group(int residueNumber, GroupInformation groupInformation, boolean parentChainIsTerminated, boolean hetAtm) {
        this.residueNumber = residueNumber;
        this.atoms = new ArrayList<>();
        this.groupInformation = groupInformation;
        this.parentChainIsTerminated = parentChainIsTerminated;
        this.hetAtm = hetAtm;
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

    public void setResidueNumber(int residueNumber) {
        this.residueNumber = residueNumber;
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
        setFeatureContainer(group.getFeatureContainer());
    }

    public Selection.AtomSelection select() {
        return Selection.on(this);
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
        return !parentChainIsTerminated && groupInformation.getType().contains("PEPTIDE LINKING");
    }

    public boolean isNucleotide() {
        return !parentChainIsTerminated && groupInformation.getType().contains("NA LINKING");
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
        return getClass().getSimpleName() + " identifier='" + getIdentifier() + "'";
    }

    /**
     * Registers a child. This object will assign a reference to itself to the atom.
     * @param atom the atom to processUniProtId
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
        return parentChainIsTerminated;
    }

    public boolean isHetAtm() {
        return hetAtm;
    }
}