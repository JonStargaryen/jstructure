package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Pyrrolysine extends AminoAcid implements NonStandardAminoAcid {
    public static final String THREE_LETTER_CODE = "PYL";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb2;
    private Atom cg2;
    private Atom cd2;
    private Atom ce2;
    private Atom n2;
    private Atom ca2;
    private Atom c2;
    private Atom o2;
    private Atom nz;
    private Atom ce;
    private Atom cd;
    private Atom cg;
    private Atom cb;

    public Pyrrolysine(Pyrrolysine pyrrolysine) {
        super(pyrrolysine);
//        this.cb2 = new Atom(pyrrolysine.cb2);
//        this.cg2 = new Atom(pyrrolysine.cg2);
//        this.cd2 = new Atom(pyrrolysine.cd2);
//        this.ce2 = new Atom(pyrrolysine.ce2);
//        this.n2 = new Atom(pyrrolysine.n2);
//        this.ca2 = new Atom(pyrrolysine.ca2);
//        this.c2 = new Atom(pyrrolysine.c2);
//        this.o2 = new Atom(pyrrolysine.o2);
//        this.nz = new Atom(pyrrolysine.nz);
//        this.ce = new Atom(pyrrolysine.ce);
//        this.cd = new Atom(pyrrolysine.cd);
//        this.cg = new Atom(pyrrolysine.cg);
//        this.cb = new Atom(pyrrolysine.cb);
    }

    public Pyrrolysine(ResidueNumber residueNumber, boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Pyrrolysine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb2() {
        return cb2;
    }

    public Atom getCg2() {
        return cg2;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getN2() {
        return n2;
    }

    public Atom getCa2() {
        return ca2;
    }

    public Atom getC2() {
        return c2;
    }

    public Atom getO2() {
        return o2;
    }

    public Atom getNz() {
        return nz;
    }

    public Atom getCe() {
        return ce;
    }

    public Atom getCd() {
        return cd;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCb() {
        return cb;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB2") && cb2 == null) {
            cb2 = atom;
        }
        if(atom.getName().equals("CG2") && cg2 == null) {
            cg2 = atom;
        }
        if(atom.getName().equals("CD2") && cd2 == null) {
            cd2 = atom;
        }
        if(atom.getName().equals("CE2") && ce2 == null) {
            ce2 = atom;
        }
        if(atom.getName().equals("N2") && n2 == null) {
            n2 = atom;
        }
        if(atom.getName().equals("CA2") && ca2 == null) {
            ca2 = atom;
        }
        if(atom.getName().equals("C2") && c2 == null) {
            c2 = atom;
        }
        if(atom.getName().equals("O2") && o2 == null) {
            o2 = atom;
        }
        if(atom.getName().equals("NZ") && nz == null) {
            nz = atom;
        }
        if(atom.getName().equals("CE") && ce == null) {
            ce = atom;
        }
        if(atom.getName().equals("CD") && cd == null) {
            cd = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
    }

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
