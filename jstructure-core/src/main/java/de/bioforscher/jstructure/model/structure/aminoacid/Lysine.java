package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Lysine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "LYS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom cd;
    private Atom ce;
    private Atom nz;

    Lysine(Lysine lysine, boolean deep) {
        super(lysine, deep);
    }

    public Lysine(ResidueIdentifier residueIdentifier,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Lysine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
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
