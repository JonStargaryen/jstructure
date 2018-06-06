package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Asparagine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ASN";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom od1;
    private Atom nd2;

    Asparagine(Asparagine asparagine, boolean deep) {
        super(asparagine, deep);
    }

    public Asparagine(ResidueIdentifier residueIdentifier,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Asparagine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getOd1() {
        return od1;
    }

    public Atom getNd2() {
        return nd2;
    }
}
