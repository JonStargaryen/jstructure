package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Lysine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "LYS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd;
    private Atom ce;
    private Atom nz;

    public Lysine(ResidueNumber residueNumber,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Lysine(ResidueNumber residueNumber) {
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

    public Atom getCe() {
        return ce;
    }

    public Atom getNz() {
        return nz;
    }
}
