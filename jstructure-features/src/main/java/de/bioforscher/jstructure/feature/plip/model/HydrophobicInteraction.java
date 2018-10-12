package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.stream.Collectors;
import java.util.stream.Stream;

public class HydrophobicInteraction extends AbstractInteraction {
//    private final Atom atom1;
//    private final Atom atom2;
//    private final Group group1;
//    private final Group group2;

    public HydrophobicInteraction(Atom atom1,
                                  Atom atom2,
                                  Group group1,
                                  Group group2) {
        super(Stream.of(atom1).collect(Collectors.toList()),
                Stream.of(atom2).collect(Collectors.toList()),
                group1,
                group2);
//        this.atom1 = atom1;
//        this.atom2 = atom2;
//        this.group1 = group1;
//        this.group2 = group2;
    }

//    public Atom getAtom1() {
//        return atom1;
//    }
//
//    public Atom getAtom2() {
//        return atom2;
//    }
//
//    public Group getGroup1() {
//        return group1;
//    }
//
//    public Group getGroup2() {
//        return group2;
//    }
}
