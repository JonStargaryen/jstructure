package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 6/19/17.
 */
public class BasicChainContainer extends AbstractBasicContainer implements ChainContainer {
    private final List<Chain> chains;
    private final List<Group> groups;
    private final List<Atom> atoms;

    public BasicChainContainer(List<Chain> chains, Protein origin) {
        super(origin);
        this.chains = chains;
        this.groups = chains.stream()
                .flatMap(Chain::groups)
                .collect(Collectors.toList());
        this.atoms = groups.stream()
                .flatMap(Group::atoms)
                .collect(Collectors.toList());
    }

    public BasicChainContainer(BasicChainContainer container) {
        super(container);
        this.chains = container.chains()
                .map(Chain::new)
                .collect(Collectors.toList());
        this.groups = chains()
                .flatMap(Chain::groups)
                .collect(Collectors.toList());
        this.atoms = groups()
                .flatMap(Group::atoms)
                .collect(Collectors.toList());
        //TODO caution with setting parent reference to null, needs some contract
        this.chains.forEach(chain -> chain.setParentProtein(null));
    }

    @Override
    public Selection.AtomSelection select() {
        return Selection.on(this);
    }

    @Override
    public LinearAlgebra.AtomContainerLinearAlgebra calculate() {
        return LinearAlgebra.on(this);
    }

    @Override
    public List<Group> getGroups() {
        return groups;
    }

    @Override
    public List<Atom> getAtoms() {
        return atoms;
    }

    @Override
    public List<Chain> getChains() {
        return chains;
    }
}