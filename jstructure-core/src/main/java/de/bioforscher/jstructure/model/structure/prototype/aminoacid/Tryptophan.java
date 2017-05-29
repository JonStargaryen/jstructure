package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Tryptophan extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "TRP";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd1;
    private Atom cd2;
    private Atom ne1;
    private Atom ce2;
    private Atom ce3;
    private Atom cz2;
    private Atom cz3;
    private Atom ch2;
    private Atom oxt;

    public Tryptophan(ResidueNumber residueNumber,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Tryptophan(ResidueNumber residueNumber) {
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

    public Atom getNe1() {
        return ne1;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getCe3() {
        return ce3;
    }

    public Atom getCz2() {
        return cz2;
    }

    public Atom getCz3() {
        return cz3;
    }

    public Atom getCh2() {
        return ch2;
    }

    public Atom getOxt() {
        return oxt;
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
        if(atom.getName().equals("CE2") && ce2 == null) {
            ce2 = atom;
        }
        if(atom.getName().equals("NE1") && ne1 == null) {
            ne1 = atom;
        }
        if(atom.getName().equals("CE3") && ce3 == null) {
            ce3 = atom;
        }
        if(atom.getName().equals("CZ2") && cz2 == null) {
            cz2 = atom;
        }
        if(atom.getName().equals("CZ3") && cz3 == null) {
            cz3 = atom;
        }
        if(atom.getName().equals("CH2") && ch2 == null) {
            ch2 = atom;
        }
        if(atom.getName().equals("OXT") && oxt == null) {
            oxt = atom;
        }
    }
}
