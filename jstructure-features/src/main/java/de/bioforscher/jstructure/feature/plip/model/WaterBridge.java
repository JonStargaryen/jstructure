package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.stream.Collectors;
import java.util.stream.Stream;

public class WaterBridge extends AbstractInteraction implements DirectedInteraction {
//    private final Atom donor;
//    private final Atom acceptor;
//    private final Group donorGroup;
//    private final Group acceptorGroup;

    public WaterBridge(Atom donor,
                       Atom acceptor,
                       Group donorGroup,
                       Group acceptorGroup) {
        super(Stream.of(donor).collect(Collectors.toList()),
                Stream.of(acceptor).collect(Collectors.toList()),
                donorGroup,
                acceptorGroup);
//        this.donor = donor;
//        this.acceptor = acceptor;
//        this.donorGroup = donorGroup;
//        this.acceptorGroup = acceptorGroup;
    }

//    public Atom getDonor() {
//        return donor;
//    }
//
//    public Atom getAcceptor() {
//        return acceptor;
//    }
//
//    public Group getDonorGroup() {
//        return donorGroup;
//    }
//
//    public Group getAcceptorGroup() {
//        return acceptorGroup;
//    }
}
