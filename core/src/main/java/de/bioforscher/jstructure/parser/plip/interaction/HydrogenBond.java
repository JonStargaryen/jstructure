package de.bioforscher.jstructure.parser.plip.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

/**
 * A hydrogen bond as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class HydrogenBond extends PLIPInteraction {
    private Atom donor, acceptor;
    private boolean sidechain, protIsDon;
    private double distanceHA, distanceDA, angle;

    HydrogenBond(Group group, Element describingElement) {
        super(group, describingElement);
        this.sidechain = getBooleanValueOfTag("sidechain");
        this.distanceHA = getDoubleValueOfTag("dist_h-a");
        this.distanceDA = getDoubleValueOfTag("dist_d-a");
        this.angle = getDoubleValueOfTag("don_angle");
        this.protIsDon = getBooleanValueOfTag("protisdon");
        this.donor = resolveAtom("donoridx");
        this.acceptor = resolveAtom("acceptoridx");
        this.partner2 = donor.getParentGroup().equals(partner1) ? acceptor.getParentGroup() : donor.getParentGroup();
    }

    public boolean isSidechain() {
        return sidechain;
    }

    public boolean isProtIsDon() {
        return protIsDon;
    }

    public double getDistanceHA() {
        return distanceHA;
    }

    public double getDistanceDA() {
        return distanceDA;
    }

    public double getAngle() {
        return angle;
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
