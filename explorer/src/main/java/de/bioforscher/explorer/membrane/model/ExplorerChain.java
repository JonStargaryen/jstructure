package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.stream.Collectors;

/**
 * One protein chain in the explorer.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerChain {
    private String id, pdb, rep;
    private List<ExplorerLigand> ligands;
    private boolean isRep;
    /* interactions */
    private List<ExplorerInteraction> halogenBonds;
    private List<ExplorerInteraction> hydrogenBonds;
    private List<ExplorerInteraction> hydrophobicInteractions;
    private List<ExplorerInteraction> metalComplexes;
    private List<ExplorerInteraction> piCationInteractions;
    private List<ExplorerInteraction> piStackings;
    private List<ExplorerInteraction> saltBridges;
    private List<ExplorerInteraction> waterBridges;

    public ExplorerChain() {
    }

    public ExplorerChain(Chain chain, String representativeChainId) {
        this.id = chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId();
        this.pdb = chain.composePDBRecord();
        this.rep = representativeChainId;
        this.isRep = representativeChainId.equals(this.id);

        this.ligands = chain.select()
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());
    }

    public String getId() {
        return id;
    }

    public String getPdb() {
        return pdb;
    }

    public String getRep() {
        return rep;
    }

    public boolean isRep() {
        return isRep;
    }

    public List<ExplorerLigand> getLigands() {
        return ligands;
    }
}
