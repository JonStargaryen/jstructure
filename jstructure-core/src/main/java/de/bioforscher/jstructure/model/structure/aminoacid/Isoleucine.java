package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Isoleucine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ILE";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg1;
    private Atom cg2;
    private Atom cd1;

    Isoleucine(Isoleucine isoleucine, boolean deep) {
        super(isoleucine, deep);
    }

    public Isoleucine(ResidueIdentifier residueIdentifier,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Isoleucine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg1() {
        return cg1;
    }

    public Atom getCg2() {
        return cg2;
    }

    public Atom getCd1() {
        return cd1;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG1") && cg1 == null) {
            cg1 = atom;
        }
        if(atom.getName().equals("CG2") && cg2 == null) {
            cg2 = atom;
        }
        if(atom.getName().equals("CD1") && cd1 == null) {
            cd1 = atom;
        }
    }
}
