package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Lysine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "LYS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd;
    private Atom ce;
    private Atom nz;

    public Lysine(Lysine lysine) {
        super(lysine);
//        this.cb = new Atom(lysine.cb);
//        this.cg = new Atom(lysine.cg);
//        this.cd = new Atom(lysine.cd);
//        this.ce = new Atom(lysine.ce);
//        this.nz = new Atom(lysine.nz);
    }

    public Lysine(ResidueIdentifier residueIdentifier,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Lysine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd() {
        return cd;
    }

    public Atom getCe() {
        return ce;
    }

    public Atom getNz() {
        return nz;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("CD") && cd == null) {
            cd = atom;
        }
        if(atom.getName().equals("CE") && ce == null) {
            ce = atom;
        }
        if(atom.getName().equals("NZ") && nz == null) {
            nz = atom;
        }
    }
}
