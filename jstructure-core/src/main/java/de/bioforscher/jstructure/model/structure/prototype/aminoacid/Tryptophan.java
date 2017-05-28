package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Tryptophan extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "TRP";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd1;
    private Atom cd2;
    private Atom ne1;
    private Atom ce2;
    private Atom ce3;
    private Atom cz2;
    private Atom cz3;
    private Atom ch2;

    public Tryptophan(ResidueNumber residueNumber,
                      boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Tryptophan(ResidueNumber residueNumber) {
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

    public Atom getNe1() {
        return ne1;
    }

    public Atom getCe2() {
        return ce2;
    }

    public Atom getCe3() {
        return ce3;
    }

    public Atom getCz2() {
        return cz2;
    }

    public Atom getCz3() {
        return cz3;
    }

    public Atom getCh2() {
        return ch2;
    }
}
