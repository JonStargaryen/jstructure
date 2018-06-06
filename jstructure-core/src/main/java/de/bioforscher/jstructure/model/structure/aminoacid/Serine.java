package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Serine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "SER";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom og;

    Serine(Serine serine, boolean deep) {
        super(serine, deep);
    }

    public Serine(ResidueIdentifier residueIdentifier,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Serine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getOg() {
        return og;
    }
}
