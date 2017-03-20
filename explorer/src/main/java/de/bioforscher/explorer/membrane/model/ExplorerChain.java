package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.stream.Collectors;

/**
 * One protein chain in the explorer.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerChain {
    private String id, sequence, pdb, rep, title;
    private List<String> homologous;
    private List<ExplorerGroup> aminoAcids;
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

    public ExplorerChain(Chain chain, List<String> homologous, String representativeChainId) {
        this.id = chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId();
        this.sequence = chain.getAminoAcidSequence();
        this.pdb = chain.composePDBRecord();
        this.rep = representativeChainId;
        this.homologous = homologous;
        this.isRep = representativeChainId.equals(this.id);
        this.title = chain.getParentProtein().getTitle();

        this.aminoAcids = chain.aminoAcids()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());

        this.ligands = chain.select()
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());
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

    public String getTitle() {
        return title;
    }

    public List<ExplorerGroup> getAminoAcids() {
        return aminoAcids;
    }

    public List<ExplorerLigand> getLigands() {
        return ligands;
    }
}
