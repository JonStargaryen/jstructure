package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxyguanosine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DG";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    public Deoxyguanosine(ResidueNumber residueNumber,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Deoxyguanosine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }
}
