package de.bioforscher.jstructure.model.structure.prototype.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.prototype.GroupPrototype;

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

    public Selenomethionine(ResidueNumber residueNumber,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Selenomethionine(ResidueNumber residueNumber) {
        this(residueNumber, false);
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
