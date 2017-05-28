package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class UnknownAminoAcid extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "UNK";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;

    public UnknownAminoAcid(ResidueNumber residueNumber,
                            boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public UnknownAminoAcid(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }
}
