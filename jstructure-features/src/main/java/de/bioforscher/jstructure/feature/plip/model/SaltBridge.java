package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SaltBridge extends AbstractInteraction {
//    private final Atom atom1;
//    private final List<Atom> saltAtoms;
//    private final Group group1;
//    private final Group saltGroup;

    public SaltBridge(Atom atom1,
                      List<Atom> saltAtoms,
                      Group group1,
                      Group saltGroup) {
        super(Stream.of(atom1).collect(Collectors.toList()),
                saltAtoms,
                group1,
                saltGroup);
//        this.atom1 = atom1;
//        this.saltAtoms = saltAtoms;
//        this.group1 = group1;
//        this.saltGroup = saltGroup;
    }

//    public Atom getAtom1() {
//        return atom1;
//    }
//
//    public List<Atom> getSaltAtoms() {
//        return saltAtoms;
//    }
//
//    public Group getGroup1() {
//        return group1;
//    }
//
//    public Group getSaltGroup() {
//        return saltGroup;
//    }
}
