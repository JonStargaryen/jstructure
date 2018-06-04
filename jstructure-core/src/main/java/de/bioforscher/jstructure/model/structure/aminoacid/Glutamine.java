package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Glutamine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "GLN";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
//    private Atom cb;
    private Atom cg;
    private Atom cd;
    private Atom oe1;
    private Atom ne2;

    Glutamine(Glutamine glutamine, boolean deep) {
        super(glutamine, deep);
    }

    public Glutamine(ResidueIdentifier residueIdentifier,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Glutamine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

//    public Atom getCb() {
//        return cb;
//    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd() {
        return cd;
    }

    public Atom getOe1() {
        return oe1;
    }

    public Atom getNe2() {
        return ne2;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
//        if(atom.getName().equals("CB") && cb == null) {
//            cb = atom;
//        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("CD") && cd == null) {
            cd = atom;
        }
        if(atom.getName().equals("OE1") && oe1 == null) {
            oe1 = atom;
        }
        if(atom.getName().equals("NE2") && ne2 == null) {
            ne2 = atom;
        }
    }
}
