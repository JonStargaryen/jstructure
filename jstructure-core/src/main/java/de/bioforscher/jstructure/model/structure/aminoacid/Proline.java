package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Proline extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "PRO";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom cd;

    Proline(Proline proline, boolean deep) {
        super(proline, deep);
    }

    public Proline(ResidueIdentifier residueIdentifier,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Proline(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd() {
        return cd;
    }
}
