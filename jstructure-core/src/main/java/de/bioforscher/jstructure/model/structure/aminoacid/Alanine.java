package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Alanine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ALA";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    Alanine(Alanine alanine, boolean deep) {
        super(alanine, deep);
    }

    public Alanine(ResidueIdentifier residueIdentifier,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Alanine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }
}
