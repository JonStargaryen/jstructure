package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * The abstract representation of an amino acid.
 * Created by bittrich on 5/24/17.
 */
public abstract class AminoAcid extends Group implements StandardAminoAcidIndicator {
    private Atom n;
    private Atom ca;
    private Atom c;
    private Atom o;
    private Atom h;

    AminoAcid(String threeLetterCode,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(threeLetterCode,
                residueNumber,
                ligand);
    }

    AminoAcid(String threeLetterCode,
              ResidueNumber residueNumber) {
        this(threeLetterCode, residueNumber, false);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(groupPrototype,
                residueNumber,
                ligand);
    }

    AminoAcid(GroupPrototype groupPrototype,
              ResidueNumber residueNumber) {
        this(groupPrototype, residueNumber, false);
    }

    public Atom getN() {
        return n;
    }

    public Atom getCa() {
        return ca;
    }

    public Atom getC() {
        return c;
    }

    public Atom getO() {
        return o;
    }

    public Atom getH() {
        return h;
    }

    public String getOneLetterCode() {
        return getGroupPrototype().getOneLetterCode().orElse("?");
    }

    @Override
    protected void addAtomInternal(Atom atom) {
        if(atom.getName().equals("N") && n == null) {
            n = atom;
        }
        if(atom.getName().equals("CA") && ca == null) {
            ca = atom;
        }
        if(atom.getName().equals("C") && c == null) {
            c = atom;
        }
        if(atom.getName().equals("O") && o == null) {
            o = atom;
        }
        if(atom.getName().equals("H") && h == null) {
            h = atom;
        }
        addSideChainAtom(atom);
    }

    protected abstract void addSideChainAtom(Atom atom);
}
