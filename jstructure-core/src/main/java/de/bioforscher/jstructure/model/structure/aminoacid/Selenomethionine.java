package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Selenomethionine extends AminoAcid implements NonStandardAminoAcid {
    public static final String THREE_LETTER_CODE = "MSE";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom se;
    private Atom ce;

    Selenomethionine(Selenomethionine selenomethionine, boolean deep) {
        super(selenomethionine, deep);
    }

    public Selenomethionine(ResidueIdentifier residueIdentifier,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Selenomethionine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    @Override
    public Class<Methionine> getParentAminoAcid() {
        return Methionine.class;
    }
}
