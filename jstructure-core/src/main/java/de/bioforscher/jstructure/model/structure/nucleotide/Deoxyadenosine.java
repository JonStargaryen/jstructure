package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxyadenosine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DA";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom n9;
    private Atom c8;
    private Atom n7;
    private Atom n6;

    public Deoxyadenosine(Deoxyadenosine deoxyadenosine) {
        super(deoxyadenosine);
    }

    public Deoxyadenosine(ResidueNumber residueNumber,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Deoxyadenosine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    @Override
    protected void addBaseAtom(Atom atom) {
        if(atom.getName().equals("N9") && n9 == null) {
            n9 = atom;
        }
        if(atom.getName().equals("C8") && c8 == null) {
            c8 = atom;
        }
        if(atom.getName().equals("N7") && n7 == null) {
            n7 = atom;
        }
        if(atom.getName().equals("N6") && n6 == null) {
            n6 = atom;
        }
    }

    public Atom getN9() {
        return n9;
    }

    public Atom getC8() {
        return c8;
    }

    public Atom getN7() {
        return n7;
    }

    public Atom getN6() {
        return n6;
    }
}
