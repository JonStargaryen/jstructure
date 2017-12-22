package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * The container element representing a {@link Structure} chain. Is composed of {@link Group} objects.
 * Created by S on 27.09.2016.
 */
public class Chain extends AbstractFeatureable implements GroupContainer {
    /**
     * reference to an undefined chain - this is used by groups without explicit parent reference
     */
    public static final Chain UNKNOWN_CHAIN = new Chain(ChainIdentifier.UNKNOWN_CHAIN_IDENTIFIER);

    private List<Group> groups;
    /**
     * The unique chain name. Usually one character, e.g. 'A'.
     */
    private ChainIdentifier chainIdentifier;
    /**
     * Handle to the containing element.
     */
    private Structure parentStructure;
    private String identifier;

    /**
     * Constructor for chain objects.
     * @param chainIdentifier the unique name of this chain
     */
    public Chain(ChainIdentifier chainIdentifier) {
        this.chainIdentifier = chainIdentifier;
        this.groups = new ArrayList<>();
    }

    Chain(Chain chain, boolean deep) {
        this.chainIdentifier = chain.chainIdentifier;
        this.identifier = chain.identifier;
        if(deep) {
            this.groups = chain.groups()
                    .map(Group::createDeepCopy)
                    .collect(Collectors.toList());
            this.groups.forEach(group -> group.setParentChain(this));
            this.parentStructure = chain.parentStructure;
        } else {
            this.groups = new ArrayList<>();
        }
    }

    public void setChainIdentifier(ChainIdentifier chainIdentifier) {
        this.chainIdentifier = chainIdentifier;
    }

    public Selection.GroupSelection select() {
        return Selection.on(this);
    }

    public LinearAlgebra.AtomContainerLinearAlgebra calculate() {
        return LinearAlgebra.on(this);
    }

    /**
     * Registers a child. This object will assign a reference to itself to the getResidue.
     * @param group the residue to process
     */
    public void addGroup(Group group) {
        getGroups().add(group);
        group.setParentChain(this);
    }

    /**
     * Returns the unique name of this chain.
     * @return a {@link ChainIdentifier}
     */
    public ChainIdentifier getChainIdentifier() {
        return chainIdentifier;
    }

    /**
     * Package-private method to set the parent reference.
     * @param parentStructure the parent
     */
    void setParentStructure(Structure parentStructure) {
        this.parentStructure = parentStructure;
    }

    /**
     * Returns the {@link Structure} this chain is associated to. If none was set, this chain points to the
     * {@link Structure#UNKNOWN_STRUCTURE}, but the unknown structure has no knowledge of the existence of this object.
     * @return the parent container
     */
    public Structure getParentStructure() {
        return parentStructure != null ? parentStructure : Structure.UNKNOWN_STRUCTURE;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " '" + getIdentifier() + "' groups='" + getGroups().size() + "'";
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? chainIdentifier.getFullName() : identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
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
    public Chain createDeepCopy() {
        return new Chain(this, true);
    }

    @Override
    public Chain createShallowCopy() {
        return new Chain(this, false);
    }
}