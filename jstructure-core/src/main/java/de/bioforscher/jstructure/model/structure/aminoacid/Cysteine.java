package de.bioforscher.jstructure.model.structure.aminoacid;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.GroupPrototype;

/**
 * Created by bittrich on 5/24/17.
 */
public class Cysteine extends AminoAcid implements StandardAminoAcid {
    public static final String THREE_LETTER_CODE = "CYS";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom cb;
    private Atom sg;

    public Cysteine(ResidueNumber residueNumber,
                    boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Cysteine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    @Override
    protected void addSideChainAtom(Atom atom) {
        if(atom.getName().equals("CB") && cb == null) {
            cb = atom;
        }
        if(atom.getName().equals("SG") && sg == null) {
            sg = atom;
        }
    }
}
