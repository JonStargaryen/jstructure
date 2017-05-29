package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

import java.util.ArrayList;
import java.util.List;

/**
 * The abstract representation of a group within a molecular structure. Some group implementations provide explicit
 * functions to retrieve specific atom. They are registered automatically during parsing. The atom registered first will
 * always be kept even when e.g. multiple alpha carbons would be provided for an amino acid. Also, these methods may
 * return null, when the corresponding atom was not provided by a PDB file. These methods are used to speed up access
 * to atoms as a more convenient, but slower alternative use the {@link de.bioforscher.jstructure.model.structure.selection.Selection}
 * implementation accessible by {@link Group#select()}
 * Created by bittrich on 5/24/17.
 */
public class Group extends AbstractFeatureable {
    private String threeLetterCode;
    private ResidueNumber residueNumber;
    private GroupPrototype groupPrototype;
    private boolean ligand;
    private List<Atom> atoms;
    private Chain parentChain;

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

    public String getThreeLetterCode() {
        return threeLetterCode;
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

    public Chain getParentChain() {
        return parentChain;
    }

    public GroupPrototype.PolymerType getPolymerType() {
        return getGroupPrototype().getPolymerType();
    }

    public boolean isAminoAcid() {
        return !isLigand() && getPolymerType() == GroupPrototype.PolymerType.PEPTIDE_LINKING;
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
//        atom.setParentGroup(this);
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
}
