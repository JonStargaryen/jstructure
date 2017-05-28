package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Threonine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "THR";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom og1;
    private Atom cg2;

    public Threonine(ResidueNumber residueNumber,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Threonine(ResidueNumber residueNumber) {
        this(residueNumber, false);
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
}
