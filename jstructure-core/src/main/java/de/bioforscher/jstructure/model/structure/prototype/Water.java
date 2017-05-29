package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Represents water molecules.
 * Created by bittrich on 5/24/17.
 */
public class Water extends Group {
    public static final String THREE_LETTER_CODE = "HOH";
    public static final GroupPrototype GROUP_PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);
    private Atom o;

    public Water(ResidueNumber residueNumber) {
        super(GROUP_PROTOTYPE,
                residueNumber,
                true);
    }

    public Atom getO() {
        return o;
    }

    @Override
    protected void addAtomInternal(Atom atom) {
        if(atom.getName().equals("O") && o == null) {
            o = atom;
        }
    }
}
