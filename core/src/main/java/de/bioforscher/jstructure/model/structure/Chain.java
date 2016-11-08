package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.GroupContainer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * The container element representing a {@link Protein} getChain. Is composed of {@link Residue} objects.
 * Created by S on 27.09.2016.
 */
public class Chain implements GroupContainer, AtomRecordWriter {
    /**
     * The container of all residues associated to this getChain.
     */
    private List<Group> groups;
    /**
     * The unique getChain name. Usually one character, e.g. 'A'.
     */
    private String chainId;
    /**
     * Handle to the containing element.
     */
    private Protein parentProtein;
    private Map<Enum, Object> featureMap;

    /**
     * Constructor for getChain objects.
     * @param chainId the unique name of this getChain
     */
    public Chain(String chainId) {
        this.chainId = chainId;
        this.groups = new ArrayList<>();
        this.featureMap = new HashMap<>();
    }

    public List<Group> getGroups() {
        return groups;
    }

    public List<Atom> getAtoms() {
        return groups().flatMap(Group::atoms)
                       .collect(Collectors.toList());
    }

    /**
     * Registers a child. This object will assign a reference to itself to the getResidue.
     * @param group the getResidue to process
     */
    public void addGroup(Group group) {
        groups.add(group);
        group.setParentChain(this);
    }

    /**
     * Returns the unique name of this getChain.
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
     * Returns the {@link Protein} this getChain is associated to.
     * @return the parent container
     */
    public Protein getParentProtein() {
        return parentProtein;
    }

    @Override
    public String toString() {
        return this.getClass().getSimpleName() + " name='" + this.chainId + "'";
    }

    @Override
    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}