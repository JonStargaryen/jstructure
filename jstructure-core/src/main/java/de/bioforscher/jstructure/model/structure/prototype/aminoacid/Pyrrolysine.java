package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Pyrrolysine extends AminoAcid implements NonStandardAminoAcid {
    public static final String THREE_LETTER_CODE = "PYL";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb2;
    private Atom cg2;
    private Atom cd2;
    private Atom ce2;
    private Atom n2;
    private Atom c2;
    private Atom o2;
    private Atom nz;
    private Atom ce;
    private Atom cd;
    private Atom cb;

    public Pyrrolysine(ResidueNumber residueNumber, boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Pyrrolysine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb2() {
        return cb2;
    }

    public Atom getCg2() {
        return cg2;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getN2() {
        return n2;
    }

    public Atom getC2() {
        return c2;
    }

    public Atom getO2() {
        return o2;
    }

    public Atom getNz() {
        return nz;
    }

    public Atom getCe() {
        return ce;
    }

    public Atom getCd() {
        return cd;
    }

    public Atom getCb() {
        return cb;
    }

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
