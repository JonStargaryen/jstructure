package de.bioforscher.jstructure.model.structure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The container element representing a {@link Protein} chain. Is composed of {@link Residue} objects.
 * Created by S on 27.09.2016.
 */
public class Chain implements GroupContainer, AtomRecordWriter {
    /**
     * The container of all residues associated to this chain.
     */
    private List<Group> groups;
    /**
     * The unique chain name. Usually one character, e.g. 'A'.
     */
    private String chainId;
    /**
     * Handle to the containing element.
     */
    private Protein parentProtein;
    private Map<String, Object> featureMap;

    /**
     * Constructor for chain objects.
     * @param chainId the unique name of this chain
     */
    public Chain(String chainId) {
        this.chainId = chainId;
        this.groups = new ArrayList<>();
        this.featureMap = new HashMap<>();
    }

    /**
     * Registers a child. This object will assign a reference to itself to the residue.
     * @param group the residue to process
     */
    public void addGroup(Group group) {
        groups.add(group);
        group.setParentChain(this);
    }

    /**
     * Returns the unique name of this chain.
     * @return a String, usually containing only a single char
     */
    public String getChainId() {
        return chainId;
    }

    /**
     * Package-private method to set the parent reference.
     * @param parentProtein the parent
     */
    void setParentProtein(Protein parentProtein) {
        this.parentProtein = parentProtein;
    }

    /**
     * Returns the {@link Protein} this chain is associated to.
     * @return the parent container
     */
    public Protein getParentProtein() {
        return parentProtein;
    }

    @Override
    public String composePDBRecord() {
        return residues().map(Residue::composePDBRecord).collect(Collectors.joining(System.lineSeparator()));
    }

    @Override
    public String toString() {
        return this.getClass().getSimpleName() + " name='" + this.chainId + "'";
    }

    public Stream<Residue> residues() {
        return groups.stream().filter(Residue.class::isInstance).map(Residue.class::cast);
    }

    @Override
    public Stream<Atom> atoms() {
        return residues().flatMap(Residue::atoms);
    }

    @Override
    public Map<String, Object> getFeatureMap() {
        return featureMap;
    }

    @Override
    public Stream<Group> groups() {
        return groups.stream();
    }
}