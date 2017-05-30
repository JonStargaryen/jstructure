package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxyadenosine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DA";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    public Deoxyadenosine(ResidueNumber residueNumber,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Deoxyadenosine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }
}
