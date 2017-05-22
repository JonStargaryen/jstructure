package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * One part of a protein around a peculiar amino acid.
 * Created by bittrich on 4/21/17.
 */
public class ExplorerMutationEnvironment {
    private String id, pdb;
    private List<ExplorerGroup> groups;
    /* interactions */
    private List<ExplorerInteraction> halogenBonds;
    private List<ExplorerInteraction> hydrogenBonds;
    private List<ExplorerInteraction> metalComplexes;
    private List<ExplorerInteraction> piCationInteractions;
    private List<ExplorerInteraction> piStackings;
    private List<ExplorerInteraction> saltBridges;
    private List<ExplorerInteraction> waterBridges;
    private Set<Integer> groupIndices;

    public ExplorerMutationEnvironment() {
    }

    public ExplorerMutationEnvironment(Chain environment, ExplorerChain explorerChain, String id) {
        this.pdb = environment.getPdbRepresentation();
        this.id = id;
        this.groups = environment.groups()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());

        this.groupIndices = groups.stream()
                .map(ExplorerGroup::getResn)
                .collect(Collectors.toSet());
        this.halogenBonds = filterRelevantInteractions(explorerChain.getHalogenBonds());
        this.hydrogenBonds = filterRelevantInteractions(explorerChain.getHydrogenBonds());
        this.metalComplexes = filterRelevantInteractions(explorerChain.getMetalComplexes());
        this.piCationInteractions = filterRelevantInteractions(explorerChain.getPiCationInteractions());
        this.piStackings = filterRelevantInteractions(explorerChain.getPiStackings());
        this.saltBridges = filterRelevantInteractions(explorerChain.getSaltBridges());
        this.waterBridges = filterRelevantInteractions(explorerChain.getWaterBridges());
    }

    private List<ExplorerInteraction> filterRelevantInteractions(List<ExplorerInteraction> interactions) {
        return interactions.stream()
                .filter(this::isRelevantInteraction)
                .collect(Collectors.toList());
    }

    private boolean isRelevantInteraction(ExplorerInteraction explorerInteraction) {
        return groupIndices.contains(explorerInteraction.getPartner1()) && groupIndices.contains(explorerInteraction.getPartner2());
    }

    public String getId() {
        return id;
    }

    public String getPdb() {
        return pdb;
    }

    public List<ExplorerGroup> getGroups() {
        return groups;
    }

    public List<ExplorerInteraction> getHalogenBonds() {
        return halogenBonds;
    }

    public List<ExplorerInteraction> getHydrogenBonds() {
        return hydrogenBonds;
    }

    public List<ExplorerInteraction> getMetalComplexes() {
        return metalComplexes;
    }

    public List<ExplorerInteraction> getPiCationInteractions() {
        return piCationInteractions;
    }

    public List<ExplorerInteraction> getPiStackings() {
        return piStackings;
    }

    public List<ExplorerInteraction> getSaltBridges() {
        return saltBridges;
    }

    public List<ExplorerInteraction> getWaterBridges() {
        return waterBridges;
    }
}
