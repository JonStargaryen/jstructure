package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Alanine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ALA";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;

    public Alanine(Alanine alanine) {
        super(alanine);
//        this.cb = new Atom(alanine.cb);
    }

    public Alanine(ResidueIdentifier residueIdentifier,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Alanine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCb() {
        return cb;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
    }
}
