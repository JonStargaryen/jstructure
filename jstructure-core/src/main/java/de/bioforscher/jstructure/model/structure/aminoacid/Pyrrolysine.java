package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

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
    private Atom ca2;
    private Atom c2;
    private Atom o2;
    private Atom nz;
    private Atom ce;
    private Atom cd;
    private Atom cg;

    Pyrrolysine(Pyrrolysine pyrrolysine, boolean deep) {
        super(pyrrolysine, deep);
    }

    public Pyrrolysine(ResidueIdentifier residueIdentifier, boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Pyrrolysine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
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

    public Atom getCa2() {
        return ca2;
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

    public Atom getCg() {
        return cg;
    }

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
