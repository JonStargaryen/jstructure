package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class Cytidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "C";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2prime;
    private Atom o2;
    private Atom n4;

    Cytidine(Cytidine cytidine, boolean deep) {
        super(cytidine, deep);
    }

    public Cytidine(ResidueIdentifier residueIdentifier,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Cytidine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getO2prime() {
        return o2prime;
    }

    public Atom getO2() {
        return o2;
    }

    public Atom getN4() {
        return n4;
    }
}
