package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Glycine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "GLY";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    Glycine(Glycine glycine, boolean deep) {
        super(glycine, deep);
    }

    public Glycine(ResidueIdentifier residueIdentifier,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Glycine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }
}
