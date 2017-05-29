package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Tyrosine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "TYR";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd1;
    private Atom cd2;
    private Atom ce1;
    private Atom ce2;
    private Atom cz;
    private Atom oh;

    public Tyrosine(ResidueNumber residueNumber,
                    boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Tyrosine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd1() {
        return cd1;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe1() {
        return ce1;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getCz() {
        return cz;
    }

    public Atom getOh() {
        return oh;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("CD1") && cd1 == null) {
            cd1 = atom;
        }
        if(atom.getName().equals("CD2") && cd2 == null) {
            cd2 = atom;
        }
        if(atom.getName().equals("CE1") && ce1 == null) {
            ce1 = atom;
        }
        if(atom.getName().equals("CE2") && ce2 == null) {
            ce2 = atom;
        }
        if(atom.getName().equals("CZ") && cz == null) {
            cz = atom;
        }
        if(atom.getName().equals("OH") && oh == null) {
            oh = atom;
        }
    }
}
