package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class AsparticAcid extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ASP";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cg;
    private Atom od1;
    private Atom od2;

    AsparticAcid(AsparticAcid asparticAcid, boolean deep) {
        super(asparticAcid, deep);
    }

    public AsparticAcid(ResidueIdentifier residueIdentifier,
                        boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public AsparticAcid(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getOd1() {
        return od1;
    }

    public Atom getOd2() {
        return od2;
    }
}
