package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxycytidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DC";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    public Deoxycytidine(ResidueNumber residueNumber,
                         boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Deoxycytidine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }
}
