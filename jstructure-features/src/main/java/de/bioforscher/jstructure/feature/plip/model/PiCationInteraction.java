package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class PiCationInteraction extends AbstractInteraction implements DirectedInteraction {
//    private final Atom cation;
//    private final List<Atom> piAtoms;
//    private final Group cationGroup;
//    private final Group piGroup;

    public PiCationInteraction(Atom cation,
                               List<Atom> piAtoms,
                               Group cationGroup,
                               Group piGroup) {
        super(Stream.of(cation).collect(Collectors.toList()),
                piAtoms,
                cationGroup,
                piGroup);
//        this.cation = cation;
//        this.piAtoms = piAtoms;
//        this.cationGroup = cationGroup;
//        this.piGroup = piGroup;
    }

//    public Atom getCation() {
//        return cation;
//    }
//
//    public List<Atom> getPiAtoms() {
//        return piAtoms;
//    }
//
//    public Group getCationGroup() {
//        return cationGroup;
//    }
//
//    public Group getPiGroup() {
//        return piGroup;
//    }
}
