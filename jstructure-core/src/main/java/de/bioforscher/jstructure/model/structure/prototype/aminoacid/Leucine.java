package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Leucine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "LEU";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd1;
    private Atom cd2;

    public Leucine(ResidueNumber residueNumber,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Leucine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd1() {
        return cd1;
    }

    public Atom getCd2() {
        return cd2;
    }
}