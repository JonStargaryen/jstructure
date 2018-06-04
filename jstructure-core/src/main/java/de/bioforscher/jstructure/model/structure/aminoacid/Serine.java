package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Serine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "SER";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
//    private Atom cb;
    private Atom og;

    Serine(Serine serine, boolean deep) {
        super(serine, deep);
    }

    public Serine(ResidueIdentifier residueIdentifier,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Serine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

//    public Atom getCb() {
//        return cb;
//    }

    public Atom getOg() {
        return og;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
//        if(atom.getName().equals("CB") && cb == null) {
//            cb = atom;
//        }
        if(atom.getName().equals("OG") && og == null) {
            og = atom;
        }
    }
}
