package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Created by bittrich on 5/30/17.
 */
public class Deoxycytidine extends Nucleotide implements StandardNucleotide {
    public static final String THREE_LETTER_CODE = "DC";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o2;
    private Atom n4;

    public Deoxycytidine(Deoxycytidine deoxycytidine) {
        super(deoxycytidine);
    }

    public Deoxycytidine(ResidueIdentifier residueIdentifier,
                         boolean ligand) {
        super(GROUP_PROTOTYPE, residueIdentifier, ligand);
    }

    public Deoxycytidine(ResidueIdentifier residueIdentifier) {
        this(residueIdentifier, false);
    }

    @Override
    protected void addBaseAtom(Atom atom) {
        if(atom.getName().equals("O2") && o2 == null) {
            o2 = atom;
        }
        if(atom.getName().equals("N4") && n4 == null) {
            n4 = atom;
        }
    }

    public Atom getO2() {
        return o2;
    }

    public Atom getN4() {
        return n4;
    }
}
