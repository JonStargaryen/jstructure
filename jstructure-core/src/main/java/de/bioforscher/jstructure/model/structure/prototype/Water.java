package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * Represents water molecules.
 * Created by bittrich on 5/24/17.
 */
public class Water extends Ligand {
    private static final String THREE_LETTER_CODE = "HOH";
    private static final Water PROTOTYPE = createPrototypeInstance(THREE_LETTER_CODE);

    private Atom o;
    private Atom h1;
    private Atom h2;

    @Override
    public Water getPrototype() {
        return PROTOTYPE;
    }
}
