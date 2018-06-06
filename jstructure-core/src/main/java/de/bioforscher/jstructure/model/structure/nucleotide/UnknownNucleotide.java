package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/30/17.
 */
public class UnknownNucleotide extends Nucleotide implements StandardNucleotide {
    //TODO better fallback for unknown nucleotides
    public static final String THREE_LETTER_CODE = "N";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    UnknownNucleotide(UnknownNucleotide unknownNucleotide, boolean deep) {
        super(unknownNucleotide, deep);
    }

    public UnknownNucleotide(String threeLetterCode,
                             ResidueIdentifier residueIdentifier,
                             boolean ligand) {
        super(threeLetterCode, residueIdentifier, ligand);
    }

    public UnknownNucleotide(String threeLetterCode,
                            ResidueIdentifier residueIdentifier) {
        this(threeLetterCode, residueIdentifier, false);
    }
}
