package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Methionine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "MET";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom sd;
    private Atom ce;

    public Methionine(Methionine methionine) {
        super(methionine);
//        this.cb = new Atom(methionine.cb);
//        this.cg = new Atom(methionine.cg);
//        this.sd = new Atom(methionine.sd);
//        this.ce = new Atom(methionine.ce);
    }

    public Methionine(ResidueIdentifier residueIdentifier,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Methionine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getSd() {
        return sd;
    }

    public Atom getCe() {
        return ce;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("SD") && sd == null) {
            sd = atom;
        }
        if(atom.getName().equals("CE") && ce == null) {
            ce = atom;
        }
    }
}
