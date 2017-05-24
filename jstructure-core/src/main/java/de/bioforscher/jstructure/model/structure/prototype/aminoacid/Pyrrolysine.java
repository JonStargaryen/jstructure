package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * Created by bittrich on 5/24/17.
 */
public class Pyrrolysine extends AminoAcid implements NonStandardAminoAcid {
    private Atom cb2;
    private Atom cg2;
    private Atom cd2;
    private Atom ce2;
    private Atom n2;
    private Atom c2;
    private Atom o2;
    private Atom nz;
    private Atom ce;
    private Atom cd;
    private Atom cb;

    @Override
    public Class<UnknownAminoAcid> getParentAminoAcid() {
        return UnknownAminoAcid.class;
    }
}
