package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Phenylalanine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "PHE";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd1;
    private Atom cd2;
    private Atom ce1;
    private Atom ce2;
    private Atom cz;

    public Phenylalanine(Phenylalanine phenylalanine) {
        super(phenylalanine);
//        this.cb = new Atom(phenylalanine.cb);
//        this.cg = new Atom(phenylalanine.cg);
//        this.cd1 = new Atom(phenylalanine.cd1);
//        this.cd2 = new Atom(phenylalanine.cd2);
//        this.ce1 = new Atom(phenylalanine.ce1);
//        this.ce2 = new Atom(phenylalanine.ce2);
//        this.cz = new Atom(phenylalanine.cz);
    }

    public Phenylalanine(ResidueIdentifier residueIdentifier,
                         boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Phenylalanine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
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
    }
}
