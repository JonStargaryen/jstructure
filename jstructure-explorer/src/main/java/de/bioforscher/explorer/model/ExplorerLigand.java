package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Represents ligands.
 * Created by bittrich on 4/6/17.
 */
public class ExplorerLigand {
    private int resn;
    private String tlc, type, name;

    public ExplorerLigand() {
    }

    public ExplorerLigand(Group ligand) {
        this.resn = ligand.getResidueNumber().getResidueNumber();
        this.tlc = ligand.getThreeLetterCode();
        this.type = tlc.equals("HOH") ? "water" : "ligand";
        this.name = ligand.getGroupPrototype().getName().replace("\"", "");
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
}
