package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/24/17.
 */
public class Proline extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "PRO";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom cd;

    public Proline(Proline proline) {
        super(proline);
        this.cb = proline.cb;
        this.cg = proline.cg;
        this.cd = proline.cd;
    }

    public Proline(ResidueNumber residueNumber,
                   boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Proline(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    public Atom getCb() {
        return cb;
    }

    public Atom getCg() {
        return cg;
    }

    public Atom getCd() {
        return cd;
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("CD") && cd == null) {
            cd = atom;
        }
    }
}
