package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/24/17.
 */
public class Selenomethionine extends AminoAcid implements NonStandardAminoAcid {
    public static final String THREE_LETTER_CODE = "MSE";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom cg;
    private Atom se;
    private Atom ce;

    public Selenomethionine(Selenomethionine selenomethionine) {
        super(selenomethionine);
//        this.cb = new Atom(selenomethionine.cb);
//        this.cg = new Atom(selenomethionine.cg);
//        this.se = new Atom(selenomethionine.se);
//        this.ce = new Atom(selenomethionine.ce);
    }

    public Selenomethionine(ResidueIdentifier residueIdentifier,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Selenomethionine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("CG") && cg == null) {
            cg = atom;
        }
        if(atom.getName().equals("SE") && se == null) {
            se = atom;
        }
        if(atom.getName().equals("CE") && ce == null) {
            ce = atom;
        }
    }

    @Override
    public Class<Methionine> getParentAminoAcid() {
        return Methionine.class;
    }
}
