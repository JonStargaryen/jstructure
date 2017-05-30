package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Thymidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "T";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    public Thymidine(ResidueNumber residueNumber,
                  boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Thymidine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }
}
