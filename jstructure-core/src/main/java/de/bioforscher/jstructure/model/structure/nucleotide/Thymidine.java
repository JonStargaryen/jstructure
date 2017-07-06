package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Thymidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "T";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2;
    private Atom o4;
    private Atom c7;

    public Thymidine(Thymidine thymidine) {
        super(thymidine);
    }

    public Thymidine(ResidueNumber residueNumber,
                     boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Thymidine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    @Override
    protected void addBaseAtom(Atom atom) {
        if(atom.getName().equals("O2") && o2 == null) {
            o2 = atom;
        }
        if(atom.getName().equals("O4") && o4 == null) {
            o4 = atom;
        }
        if(atom.getName().equals("C7") && c7 == null) {
            c7 = atom;
        }
    }
}
