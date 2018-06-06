package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Methionine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "MET";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom sd;
    private Atom ce;

    Methionine(Methionine methionine, boolean deep) {
        super(methionine, deep);
    }

    public Methionine(ResidueIdentifier residueIdentifier,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Methionine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getSd() {
        return sd;
    }

    public Atom getCe() {
        return ce;
    }
}
