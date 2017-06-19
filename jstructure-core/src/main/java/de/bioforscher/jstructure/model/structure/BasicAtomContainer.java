package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by bittrich on 6/19/17.
 */
public class BasicAtomContainer extends AbstractBasicContainer implements AtomContainer {
    private final List<Atom> atoms;

    public BasicAtomContainer(List<Atom> atoms, Protein origin) {
        super(origin);
        this.atoms = atoms;
    }

    public BasicAtomContainer(BasicAtomContainer container) {
        super(container);
        this.atoms = container.atoms()
                .map(Atom::new)
                .collect(Collectors.toList());
        //TODO caution with setting parent reference to null, needs some contract
        this.atoms().forEach(atom -> atom.setParentGroup(null));
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
    public List<Atom> getAtoms() {
        return atoms;
    }
}
