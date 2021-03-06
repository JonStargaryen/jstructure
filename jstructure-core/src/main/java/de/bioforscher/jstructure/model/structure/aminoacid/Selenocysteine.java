package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Selenocysteine extends AminoAcid implements NonStandardAminoAcid {
    public static final String THREE_LETTER_CODE = "SEC";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom se;

    Selenocysteine(Selenocysteine selenocysteine, boolean deep) {
        super(selenocysteine, deep);
    }

    public Selenocysteine(ResidueIdentifier residueIdentifier,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Selenocysteine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getSe() {
        return se;
    }

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
