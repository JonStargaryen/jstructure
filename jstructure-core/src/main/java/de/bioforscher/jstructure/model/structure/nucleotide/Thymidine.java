package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class Thymidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "T";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2;
    private Atom o4;
    private Atom c7;

    Thymidine(Thymidine thymidine, boolean deep) {
        super(thymidine, deep);
    }

    public Thymidine(ResidueIdentifier residueIdentifier,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Thymidine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }
}
