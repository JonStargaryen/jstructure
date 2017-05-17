package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

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

    @Override
    public String toString() {
        return getClass().getSimpleName() + " donor=" + toString(donor) + " acceptor=" + toString(acceptor);
    }
}
