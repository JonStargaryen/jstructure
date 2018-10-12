package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class PiStackingInteraction extends AbstractInteraction {
//    private final Atom atom1;
//    private final List<Atom> piAtoms;
//    private final Group group1;
//    private final Group piGroup;

    public PiStackingInteraction(Atom atom1,
                                 List<Atom> piAtoms,
                                 Group group1,
                                 Group piGroup) {
        super(Stream.of(atom1).collect(Collectors.toList()),
                piAtoms,
                group1,
                piGroup);
//        this.atom1 = atom1;
//        this.piAtoms = piAtoms;
//        this.group1 = group1;
//        this.piGroup = piGroup;
    }

//    public Atom getAtom1() {
//        return atom1;
//    }
//
//    public List<Atom> getPiAtoms() {
//        return piAtoms;
//    }
//
//    public Group getGroup1() {
//        return group1;
//    }
//
//    public Group getPiGroup() {
//        return piGroup;
//    }
}
