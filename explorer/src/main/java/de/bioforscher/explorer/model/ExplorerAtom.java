package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Atom;

/**
 * The reduced representation of {@link Atom} objects.
 * Created by bittrich on 2/22/17.
 */
public class ExplorerAtom {
    private String pdb;

    public ExplorerAtom() {

    }

    public ExplorerAtom(Atom atom) {
        this.pdb = atom.composePDBRecord();
    }

    public String getPdb() {
        return pdb;
    }
}
