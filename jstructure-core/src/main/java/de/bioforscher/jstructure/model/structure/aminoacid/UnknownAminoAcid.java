package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class UnknownAminoAcid extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "UNK";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;

    UnknownAminoAcid(UnknownAminoAcid unknownAminoAcid, boolean deep) {
        super(unknownAminoAcid, deep);
    }

    public UnknownAminoAcid(String threeLetterCode,
                            ResidueIdentifier residueIdentifier,
                            boolean ligand) {
        super(threeLetterCode, residueIdentifier, ligand);
    }

    public UnknownAminoAcid(String threeLetterCode,
                            ResidueIdentifier residueIdentifier) {
        this(threeLetterCode, residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }
}
