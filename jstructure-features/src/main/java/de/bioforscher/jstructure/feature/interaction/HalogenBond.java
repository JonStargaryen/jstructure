package de.bioforscher.jstructure.feature.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

import java.util.stream.Stream;

/**
 * A halogen bond as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class HalogenBond extends PLIPInteraction {
    private double distance;
    private double donorAngle;
    private double acceptorAngle;
    private Atom donor;
    private Atom acceptor;

    HalogenBond(Group group, Element describingElement) {
        super(group, describingElement);
        this.distance = getDoubleValueOfTag("dist");
        this.donorAngle = getDoubleValueOfTag("don_angle");
        this.acceptorAngle = getDoubleValueOfTag("acc_angle");
        this.donor = resolveAtom("don_idx");
        this.acceptor = resolveAtom("acc_idx");
        this.partner2 = donor.getParentGroup().equals(partner1) ? acceptor.getParentGroup() : donor.getParentGroup();
    }

    public double getDistance() {
        return distance;
    }

    public double getDonorAngle() {
        return donorAngle;
    }

    public double getAcceptorAngle() {
        return acceptorAngle;
    }

    public Atom getDonor() {
        return donor;
    }

    public Atom getAcceptor() {
        return acceptor;
    }

    /**
     * True iff this is an interaction to another amino acid.
     * @return true when <code>partner2.isAminoAcid()</code> evaluates to <code>true</code>
     */
    public boolean partnerIsAminoAcid() {
        return partner2.isAminoAcid();
    }

    @Override
    boolean isSane() {
        return acceptor != null && donor != null;
    }

    @Override
    public Stream<Atom> allAtoms() {
        return Stream.of(acceptor, donor);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " donor=" + toString(donor) + " acceptor=" + toString(acceptor);
    }
}
