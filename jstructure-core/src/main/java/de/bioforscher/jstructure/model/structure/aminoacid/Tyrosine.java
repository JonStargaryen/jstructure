package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Tyrosine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "TYR";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom cd1;
    private Atom cd2;
    private Atom ce1;
    private Atom ce2;
    private Atom cz;
    private Atom oh;

    Tyrosine(Tyrosine tyrosine, boolean deep) {
        super(tyrosine, deep);
    }

    public Tyrosine(ResidueIdentifier residueIdentifier,
                    boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Tyrosine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd1() {
        return cd1;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe1() {
        return ce1;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getCz() {
        return cz;
    }

    public Atom getOh() {
        return oh;
    }
}
