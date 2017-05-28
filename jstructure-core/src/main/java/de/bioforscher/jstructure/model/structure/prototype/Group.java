package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

import java.util.ArrayList;
import java.util.List;

/**
 * The abstract representation of a group within a molecular structure.
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

    protected static GroupPrototype createPrototypeInstance(String id) {
        return GroupPrototypeParser.getInstance().getPrototype(id);
    }
}
