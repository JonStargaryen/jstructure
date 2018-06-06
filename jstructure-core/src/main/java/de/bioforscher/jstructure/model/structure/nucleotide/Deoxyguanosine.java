package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxyguanosine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DG";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom n9;
    private Atom c8;
    private Atom n7;
    private Atom o6;
    private Atom n2;

    Deoxyguanosine(Deoxyguanosine deoxyguanosine, boolean deep) {
        super(deoxyguanosine, deep);
    }

    public Deoxyguanosine(ResidueIdentifier residueIdentifier,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Deoxyguanosine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getN9() {
        return n9;
    }

    public Atom getC8() {
        return c8;
    }

    public Atom getN7() {
        return n7;
    }

    public Atom getO6() {
        return o6;
    }

    public Atom getN2() {
        return n2;
    }
}
