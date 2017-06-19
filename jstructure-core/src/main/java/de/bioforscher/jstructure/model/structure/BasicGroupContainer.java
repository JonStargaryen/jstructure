package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 6/19/17.
 */
public class BasicGroupContainer extends AbstractBasicContainer implements GroupContainer {
    private final List<Group> groups;
    private final List<Atom> atoms;

    public BasicGroupContainer(List<Group> groups, Protein origin) {
        super(origin);
        this.groups = groups;
        this.atoms = groups.stream()
                .flatMap(Group::atoms)
                .collect(Collectors.toList());
    }

    public BasicGroupContainer(BasicGroupContainer container) {
        super(container);
        this.groups = container.groups()
                .map(Group::new)
                .collect(Collectors.toList());
        this.atoms = groups()
                .flatMap(Group::atoms)
                .collect(Collectors.toList());
        //TODO caution with setting parent reference to null, needs some contract
        this.groups.forEach(group -> group.setParentChain(null));
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
}