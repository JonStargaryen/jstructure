package de.bioforscher.explorer.membrane.model.representative;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * The reduced representation of a ligand.
 * Created by bittrich on 2/28/17.
 */
public class ExplorerLigand {
    private int resn;
    private String chain, tlc, type, name;

    public ExplorerLigand() {

    }

    public ExplorerLigand(Group group) {
        this.resn = group.getResidueNumber();
        this.chain = group.getParentChain().getChainId();
        this.tlc = group.getThreeLetterCode();
        if(this.tlc.equals("HOH")) {
            this.type = "water";
        } else {
            this.type = "ligand";
        }
        this.name = group.getGroupInformation().getName().replace("\"", "");
    }

    public int getResn() {
        return resn;
    }

    public String getChain() {
        return chain;
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
