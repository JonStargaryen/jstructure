package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * The container element representing a {@link Protein} getChain. Is composed of {@link Group} objects.
 * Created by S on 27.09.2016.
 */
public class Chain extends AbstractFeatureContainer implements GroupContainer {
    /**
     * reference to an undefined chain - this is used by groups without explicit parent reference
     */
    static final Chain UNKNOWN_CHAIN = new Chain("X");

    private List<Group> groups;
    /**
     * The unique getChain name. Usually one character, e.g. 'A'.
     */
    private String chainId;
    /**
     * Handle to the containing element.
     */
    private Protein parentProtein;
    private String identifier;

    /**
     * Constructor for getChain objects.
     * @param chainId the unique name of this getChain
     */
    public Chain(String chainId) {
        this.chainId = chainId;
        this.groups = new ArrayList<>();
    }

    public Chain(Chain chain) {
        // deep clone entries
        this.groups = chain.groups()
                .map(Group::new)
                .collect(Collectors.toList());
        this.groups.forEach(group -> group.setParentChain(this));
        this.chainId = chain.chainId;
        // reference parent
        this.parentProtein = chain.parentProtein;
    }

    public Chain(List<Group> groups) {
        this.groups = groups;
    }

    Chain() {

    }

    public Selection.GroupSelection select() {
        return Selection.on(this);
    }

    /**
     * Registers a child. This object will assign a reference to itself to the getResidue.
     * @param group the getResidue to process
     */
    public void addGroup(Group group) {
        getGroups().add(group);
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
        return parentProtein != null ? parentProtein : Protein.UNKNOWN_PROTEIN;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " identifier='" + getIdentifier() + "' groups='" + getGroups().size() + "'";
    }

    @Override
    public List<Group> getGroups() {
        return groups;
    }

    @Override
    public List<Atom> getAtoms() {
        return groups().flatMap(Group::atoms).collect(Collectors.toList());
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? chainId : identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}