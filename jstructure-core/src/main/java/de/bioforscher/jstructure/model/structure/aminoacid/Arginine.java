package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Arginine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ARG";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom cd;
    private Atom ne;
    private Atom cz;

    Arginine(Arginine arginine, boolean deep) {
        super(arginine, deep);
    }

    public Arginine(ResidueIdentifier residueIdentifier,
                    boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Arginine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd() {
        return cd;
    }

    public Atom getNe() {
        return ne;
    }

    public Atom getCz() {
        return cz;
    }
}
