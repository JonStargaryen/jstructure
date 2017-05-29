package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Arginine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ARG";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd;
    private Atom ne;
    private Atom cz;

    public Arginine(ResidueNumber residueNumber,
                    boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Arginine(ResidueNumber residueNumber) {
        this(residueNumber, false);
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

    public Atom getNe() {
        return ne;
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
        if(atom.getName().equals("CD") && cd == null) {
            cd = atom;
        }
        if(atom.getName().equals("NE") && ne == null) {
            ne = atom;
        }
        if(atom.getName().equals("CZ") && cz == null) {
            cz = atom;
        }
    }
}
