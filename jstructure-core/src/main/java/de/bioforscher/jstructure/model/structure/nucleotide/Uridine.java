package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class Uridine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "U";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2prime;
    private Atom o2;
    private Atom o4;

    Uridine(Uridine uridine, boolean deep) {
        super(uridine, deep);
    }

    public Uridine(ResidueIdentifier residueIdentifier,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Uridine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getO2prime() {
        return o2prime;
    }

    public Atom getO2() {
        return o2;
    }

    public Atom getO4() {
        return o4;
    }
}
