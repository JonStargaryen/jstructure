package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Histidine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "HIS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom nd1;
    private Atom cd2;
    private Atom ce1;
    private Atom ne2;

    public Histidine(ResidueNumber residueNumber,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Histidine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getNd1() {
        return nd1;
    }

    public Atom getCd2() {
        return cd2;
    }

    public Atom getCe1() {
        return ce1;
    }

    public Atom getNe2() {
        return ne2;
    }
}