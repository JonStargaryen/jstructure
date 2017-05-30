package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class UnknownNucleotide extends Nucleotide implements StandardNucleotide {
    public UnknownNucleotide(String threeLetterCode,
                            ResidueNumber residueNumber,
                            boolean ligand) {
        super(threeLetterCode, residueNumber, ligand);
    }

    public UnknownNucleotide(String threeLetterCode,
                            ResidueNumber residueNumber) {
        this(threeLetterCode, residueNumber, false);
    }
}
