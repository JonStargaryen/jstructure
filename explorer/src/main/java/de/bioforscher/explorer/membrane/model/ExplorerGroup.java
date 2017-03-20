package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Represents one amino acid.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerGroup {
    private int resn;
    private String olc, tlc;

    public ExplorerGroup() {
    }

    public ExplorerGroup(Group group) {
        this.resn = group.getResidueNumber();
        this.olc = group.getGroupInformation().getOneLetterCode();
        this.tlc = group.getThreeLetterCode();
    }

    public int getResn() {
        return resn;
    }

    public String getOlc() {
        return olc;
    }

    public String getTlc() {
        return tlc;
    }
}
