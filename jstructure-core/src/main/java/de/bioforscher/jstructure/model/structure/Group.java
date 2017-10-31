package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The generic representation of a group within a molecular structure. Some group implementations provide explicit
 * functions to retrieve specific atom. They are registered automatically during parsing. The atom registered first will
 * always be kept even when e.g. multiple alpha carbons would be provided for an amino acid. Also, these methods may
 * return null, when the corresponding atom was not provided by a PDB file. These methods are used to speed up access
 * to atoms as to a more convenient, but slower alternative in the {@link Selection} implementation accessible by
 * {@link Group#select()}.
 * Created by bittrich on 5/24/17.
 */
public class Group extends AbstractFeatureable implements AtomContainer {
    /**
     * reference to an undefined group - this is used by atoms without explicit parent reference
     */
    public static final Group UNKNOWN_GROUP = new Group("UNK",
            IdentifierFactory.createResidueIdentifier(0),
            false);
    public static final Set<String> HYDROGEN_NAMES = Stream.of("H", "D", "T").collect(Collectors.toSet());
    private String threeLetterCode;
    private ResidueIdentifier residueIdentifier;
    private GroupPrototype groupPrototype;
    private boolean ligand;
    private List<Atom> atoms;
    private Chain parentChain;
    private String identifier;

    public Group(String threeLetterCode,
                 ResidueIdentifier residueIdentifier,
                 boolean ligand) {
        this(createPrototypeInstance(threeLetterCode),
                residueIdentifier,
                ligand);
        // safety-net: maybe the group prototype cannot be created, still keep given threeLetterCode
        this.threeLetterCode = threeLetterCode;
    }

    public Group(GroupPrototype groupPrototype,
                 ResidueIdentifier residueIdentifier,
                 boolean ligand) {
        this.threeLetterCode = groupPrototype.getThreeLetterCode();
        this.residueIdentifier = residueIdentifier;
        this.groupPrototype = groupPrototype;
        this.ligand = ligand;
        this.atoms = new ArrayList<>();
    }

    protected Group(Group group, boolean deep) {
        this.threeLetterCode = group.threeLetterCode;
        this.residueIdentifier = group.residueIdentifier;
        this.groupPrototype = group.groupPrototype;
        this.ligand = group.ligand;
        this.identifier = group.identifier;
        if(deep) {
            this.atoms = group.atoms()
                    .map(Atom::createDeepCopy)
                    .collect(Collectors.toList());
            this.atoms().forEach(atom -> atom.setParentGroup(this));
            this.parentChain = group.parentChain;
        } else {
            this.atoms = new ArrayList<>();
        }
    }

    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    /**
     * Package-private method to set this group's three-letter-code. Normally this is inferred from the provided group
     * prototype, but this may lead to consistency when dealing with unknown ligands whose name was actually parsed from
     * the PDB file.
     * @param threeLetterCode the three-letter-code forced on this group
     */
    void setThreeLetterCode(String threeLetterCode) {
        this.threeLetterCode = threeLetterCode;
    }

    /**
     * Returns the residue identifier of this group, i.e. the synthesis of residue number and potential insertion code.
     * This is the naming used in the PDB format and the preferred and stable way to select individual residues. For the
     * numbering in the internal data structure use {@link #getResidueIndex()}.
     * @return the object describing this residue's identifier
     */
    public ResidueIdentifier getResidueIdentifier() {
        return residueIdentifier;
    }

    /**
     * Returns the index of this group in the parent chain. This is the plain representation of the ordering in the
     * underlying data structure. For identifiers in the PDB format use {@link #getResidueIdentifier()}.
     * @return the index of this residue in the parent chain, starting with 0
     */
    public int getResidueIndex() {
        return parentChain.getGroups().indexOf(this);
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
        // delegate to internal implementation, so that concrete impls such as an AminoAcid can infer their fields/getters correctly
        addAtomInternal(atom);
    }

    /**
     * If the child class supports specific atom identified by name.
     * @param atom the atom to be handled - identified by its name, assigned to a particular field of the child class
     */
    protected void addAtomInternal(Atom atom) {

    }

    protected static GroupPrototype createPrototypeInstance(String id) {
        return GroupPrototypeParser.getInstance().getPrototype(id);
    }

    /**
     * Returns the {@link Chain} this getResidue is associated to. If none was set, this group points to the
     * {@link Chain#UNKNOWN_CHAIN}, but the unknown chain has no knowledge of the existence of this object.
     * @return the parent container
     */
    public Chain getParentChain() {
        return parentChain != null ? parentChain : Chain.UNKNOWN_CHAIN;
    }

    @Override
    public Structure getParentStructure() {
        return getParentChain().getParentStructure();
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
        return identifier == null ? threeLetterCode + "-" + residueIdentifier : identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    @Override
    public Group createDeepCopy() {
        return callCopyConstructor(true);
    }

    @Override
    public Group createShallowCopy() {
        return callCopyConstructor(false);
    }

    private Group callCopyConstructor(boolean deep) {
        try {
            Constructor<? extends Group> constructor = getClass().getDeclaredConstructor(getClass(), boolean.class);
            constructor.setAccessible(true);
            return constructor.newInstance(this, deep);
        } catch (NoSuchMethodException | IllegalAccessException | InstantiationException | InvocationTargetException e) {
            throw new UnsupportedOperationException(e);
        }
    }
}
