package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The abstract representation of a group within a molecular structure. Some group implementations provide explicit
 * functions to retrieve specific atom. They are registered automatically during parsing. The atom registered first will
 * always be kept even when e.g. multiple alpha carbons would be provided for an amino acid. Also, these methods may
 * return null, when the corresponding atom was not provided by a PDB file. These methods are used to speed up access
 * to atoms as a more convenient, but slower alternative use the {@link Selection} implementation accessible by
 * {@link Group#select()}.
 * Created by bittrich on 5/24/17.
 */
public class Group extends AbstractFeatureable implements AtomContainer {
    /**
     * reference to an undefined group - this is used by atoms without explicit parent reference
     */
    static final Group UNKNOWN_GROUP = new Group("UNK",
            new ResidueNumber(0),
            false);
    public static final Set<String> HYDROGEN_NAMES = Stream.of("H", "D", "T").collect(Collectors.toSet());
    private String threeLetterCode;
    private ResidueNumber residueNumber;
    private GroupPrototype groupPrototype;
    private boolean ligand;
    private List<Atom> atoms;
    private Chain parentChain;
    private String identifier;

    public Group(String threeLetterCode,
                 ResidueNumber residueNumber,
                 boolean ligand) {
        this(createPrototypeInstance(threeLetterCode),
                residueNumber,
                ligand);
        // safety-net: maybe the group prototype cannot be created, still keep given threeLetterCode
        this.threeLetterCode = threeLetterCode;
    }

    public Group(GroupPrototype groupPrototype,
                 ResidueNumber residueNumber,
                 boolean ligand) {
        this.threeLetterCode = groupPrototype.getThreeLetterCode();
        this.residueNumber = residueNumber;
        this.groupPrototype = groupPrototype;
        this.ligand = ligand;
        this.atoms = new ArrayList<>();
        this.parentChain = Chain.UNKNOWN_CHAIN;
    }

    public Group(Group group) {
        this.threeLetterCode = group.threeLetterCode;
        this.residueNumber = group.residueNumber;
        this.groupPrototype = group.groupPrototype;
        this.ligand = group.ligand;
        // deep clone entries
        this.atoms = group.atoms()
                .map(Atom::new)
                .collect(Collectors.toList());
        this.atoms().forEach(atom -> atom.setParentGroup(this));
        // reference parent
        this.parentChain = group.parentChain;
        this.identifier = group.identifier;
        // set reference to feature map
        setFeatureContainer(group.getFeatureContainer());
    }

    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    public void setThreeLetterCode(String threeLetterCode) {
        this.threeLetterCode = threeLetterCode;
    }

    public ResidueNumber getResidueNumber() {
        return residueNumber;
    }

    public GroupPrototype getGroupPrototype() {
        return groupPrototype;
    }

    public List<Atom> getAtoms() {
        return atoms;
    }

    public GroupPrototype.PolymerType getPolymerType() {
        return getGroupPrototype().getPolymerType();
    }

    public boolean isAminoAcid() {
        return !isLigand() && (getPolymerType() == GroupPrototype.PolymerType.PEPTIDE_LINKING || getPolymerType() ==
                GroupPrototype.PolymerType.PEPTIDE_LIKE || getPolymerType() ==
                GroupPrototype.PolymerType.PEPTIDE_TERMINUS);
    }

    public boolean isNucleotide() {
        return !isLigand() && getPolymerType() == GroupPrototype.PolymerType.NA_LINKING;
    }

    public boolean isWater() {
        return getThreeLetterCode().equals(Water.THREE_LETTER_CODE);
    }

    public boolean isLigand() {
        return ligand;
    }

    public void addAtom(Atom atom) {
        atoms.add(atom);
        // set reference to this as parent
        atom.setParentGroup(this);
        // delegate to internal implementation
        addAtomInternal(atom);
    }

    /**
     * If the child class supports specific atom identified by name.
     * @param atom the atom to be handled - identified by its name, assigned to a particular field of the child class
     */
    protected void addAtomInternal(Atom atom) {
        //TODO this could be realized by reflection
    }

    protected static GroupPrototype createPrototypeInstance(String id) {
        return GroupPrototypeParser.getInstance().getPrototype(id);
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

    public Selection.AtomSelection select() {
        return Selection.on(this);
    }

    public LinearAlgebra.AtomContainerLinearAlgebra calculate() {
        return LinearAlgebra.on(this);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " '" + getIdentifier() + "'";
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? threeLetterCode + "-" + residueNumber : identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}
