package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Represents ligands.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerLigand {
    private int resn;
    private String tlc, type, name, pdb;

    public ExplorerLigand() {
    }

    public ExplorerLigand(Group ligand) {
        this.resn = ligand.getResidueNumber();
        this.tlc = ligand.getThreeLetterCode();
        this.type = tlc.equals("HOH") ? "water" : "ligand";
        this.name = ligand.getGroupInformation().getName().replace("\"", "");
        this.pdb = ligand.composePDBRecord();
    }

    public int getResn() {
        return resn;
    }

    public String getTlc() {
        return tlc;
    }

    public String getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    public String getPdb() {
        return pdb;
    }
}
