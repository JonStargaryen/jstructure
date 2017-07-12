package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class UnknownNucleotide extends Nucleotide implements StandardNucleotide {
    //TODO better fallback for unknown nucleotides
    public static final String THREE_LETTER_CODE = "N";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

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
