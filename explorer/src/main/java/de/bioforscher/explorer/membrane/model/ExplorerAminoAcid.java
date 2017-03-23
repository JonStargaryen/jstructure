package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Sequence positions in the msa.
 * Created by bittrich on 3/21/17.
 */
public class ExplorerAminoAcid {
    private int resn;
    private String olc;

    public ExplorerAminoAcid() {
    }

    public ExplorerAminoAcid(int resn, String olc) {
        this.resn = resn;
        this.olc = olc;
    }

    public ExplorerAminoAcid(int resn, String olc, Group aminoAcid) {
        this(resn, olc);
        //TODO assign uniprot annotations here
    }

    public int getResn() {
        return resn;
    }

    public String getOlc() {
        return olc;
    }
}
