package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class AsparticAcid extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "ASP";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom od1;
    private Atom od2;

    public AsparticAcid(ResidueNumber residueNumber,
                        boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public AsparticAcid(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
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