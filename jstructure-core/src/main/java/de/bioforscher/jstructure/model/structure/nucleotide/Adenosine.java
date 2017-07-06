package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public class Adenosine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "A";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2prime;
    private Atom n9;
    private Atom c8;
    private Atom n7;
    private Atom n6;

    public Adenosine(Adenosine adenosine) {
        super(adenosine);
    }

    public Adenosine(ResidueNumber residueNumber,
                          boolean ligand) {
        super(GROUP_PROTOTYPE, residueNumber, ligand);
    }

    public Adenosine(ResidueNumber residueNumber) {
        this(residueNumber, false);
    }

    @Override
    protected void addBaseAtom(Atom atom) {
        if(atom.getName().equals("\"O2'\"") && o2prime == null) {
            o2prime = atom;
        }
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

    public Atom getO2prime() {
        return o2prime;
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
