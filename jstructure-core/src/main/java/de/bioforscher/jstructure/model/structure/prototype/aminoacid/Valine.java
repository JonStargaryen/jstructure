package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Valine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "VAL";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg1;
    private Atom cg2;

    public Valine(ResidueNumber residueNumber,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Valine(ResidueNumber residueNumber) {
        this(residueNumber, false);
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
}
