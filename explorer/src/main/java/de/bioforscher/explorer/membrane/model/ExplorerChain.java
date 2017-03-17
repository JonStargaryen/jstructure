package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;

/**
 * One protein chain in the explorer.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerChain {
    private String id, sequence, pdb, rep;
    private List<String> homologous;
    private boolean isRep;

    public ExplorerChain() {
    }

    public ExplorerChain(Chain chain, List<String> homologous, String representativeChainId) {
        this.id = chain.getParentProtein().getName().toUpperCase() + "." + chain.getChainId();
        this.sequence = chain.getAminoAcidSequence();
        this.pdb = chain.composePDBRecord();
        this.rep = representativeChainId;
        this.homologous = homologous;
        this.isRep = representativeChainId.equals(this.id);
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }

    public String getPdb() {
        return pdb;
    }

    public List<String> getHomologous() {
        return homologous;
    }

    public String getRep() {
        return rep;
    }

    public boolean isRep() {
        return isRep;
    }
}
