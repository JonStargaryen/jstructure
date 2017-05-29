package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class GlutamicAcid extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "GLU";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd;
    private Atom oe1;
    private Atom oe2;

    public GlutamicAcid(ResidueNumber residueNumber,
                        boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public GlutamicAcid(ResidueNumber residueNumber) {
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

    public Atom getOe1() {
        return oe1;
    }

    public Atom getOe2() {
        return oe2;
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
        if(atom.getName().equals("OE1") && oe1 == null) {
            oe1 = atom;
        }
        if(atom.getName().equals("OE2") && oe2 == null) {
            oe2 = atom;
        }
    }
}
