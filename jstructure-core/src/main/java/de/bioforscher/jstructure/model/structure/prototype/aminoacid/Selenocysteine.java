package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * Created by bittrich on 5/24/17.
 */
public class Selenocysteine extends AminoAcid implements NonStandardAminoAcid {
    private Atom cb;
    private Atom se;

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
