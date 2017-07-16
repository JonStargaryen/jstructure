package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Threonine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "THR";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom og1;
    private Atom cg2;

    Threonine(Threonine threonine, boolean deep) {
        super(threonine, deep);
    }

    public Threonine(ResidueIdentifier residueIdentifier,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Threonine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getOg1() {
        return og1;
    }

    public Atom getCg2() {
        return cg2;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("OG1") && og1 == null) {
            og1 = atom;
        }
        if(atom.getName().equals("CG2") && cg2 == null) {
            cg2 = atom;
        }
    }
}
