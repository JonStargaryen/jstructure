package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Valine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "VAL";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg1;
    private Atom cg2;

    Valine(Valine valine, boolean deep) {
        super(valine, deep);
    }

    public Valine(ResidueIdentifier residueIdentifier,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Valine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg1() {
        return cg1;
    }

    public Atom getCg2() {
        return cg2;
    }
}
