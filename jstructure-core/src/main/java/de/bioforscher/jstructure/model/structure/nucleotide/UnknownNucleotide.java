package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class UnknownNucleotide extends Nucleotide implements StandardNucleotide {
    public UnknownNucleotide(String threeLetterCode,
                            ResidueIdentifier residueIdentifier,
                            boolean ligand) {
        super(threeLetterCode, residueIdentifier, ligand);
    }

    public UnknownNucleotide(String threeLetterCode,
                            ResidueIdentifier residueIdentifier) {
        this(threeLetterCode, residueIdentifier, false);
    }

    @Override
    protected void addBaseAtom(Atom atom) {

    }
}
