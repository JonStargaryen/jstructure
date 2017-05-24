package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * Created by bittrich on 5/24/17.
 */
public class Selenomethionine extends AminoAcid implements NonStandardAminoAcid {
    private Atom cb;
    private Atom cg;
    private Atom se;
    private Atom ce;

    @Override
    public Class<Methionine> getParentAminoAcid() {
        return Methionine.class;
    }
}
