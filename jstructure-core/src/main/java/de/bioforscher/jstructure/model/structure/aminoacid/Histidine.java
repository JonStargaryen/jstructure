package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Histidine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "HIS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
//    private Atom cb;
    private Atom cg;
    private Atom nd1;
    private Atom cd2;
    private Atom ce1;
    private Atom ne2;

    Histidine(Histidine histidine, boolean deep) {
        super(histidine, deep);
    }

    public Histidine(ResidueIdentifier residueIdentifier,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Histidine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

//    public Atom getCb() {
//        return cb;
//    }

    public Atom getCg() {
        return cg;
    }

    public Atom getNd1() {
        return nd1;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe1() {
        return ce1;
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
        if(atom.getName().equals("ND1") && nd1 == null) {
            nd1 = atom;
        }
        if(atom.getName().equals("CD2") && cd2 == null) {
            cd2 = atom;
        }
        if(atom.getName().equals("CE1") && ce1 == null) {
            ce1 = atom;
        }
        if(atom.getName().equals("NE2") && ne2 == null) {
            ne2 = atom;
        }
    }
}
