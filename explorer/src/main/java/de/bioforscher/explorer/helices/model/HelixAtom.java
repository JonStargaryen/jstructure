package de.bioforscher.explorer.helices.model;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * Created by bittrich on 3/7/17.
 */
public class HelixAtom {
    private String pdb;

    public HelixAtom() {

    }

    public HelixAtom(Atom atom) {
        this.pdb = atom.composePDBRecord();
    }

    public String getPdb() {
        return pdb;
    }
}
