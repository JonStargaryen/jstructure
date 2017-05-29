package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Glycine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "GLY";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    public Glycine(ResidueNumber residueNumber,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Glycine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    @Override
    protected void addSideChainAtom(Atom atom) {

    }
}
