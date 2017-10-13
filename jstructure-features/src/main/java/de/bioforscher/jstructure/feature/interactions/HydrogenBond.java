package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.nodes.Element;

import java.util.stream.Stream;

/**
 * A hydrogen bond as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class HydrogenBond extends PLIPInteraction {
    private Atom donor;
    private Atom acceptor;
    private boolean sidechain;
    private boolean protIsDon;
    private double distanceHA;
    private double distanceDA;
    private double angle;

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
    public boolean isBackboneInteraction() {
        return AminoAcid.isBackboneAtom(acceptor) && AminoAcid.isBackboneAtom(donor);
    }

    @Override
    public boolean isSideChainInteraction() {
        return AminoAcid.isSideChainAtom(acceptor) && AminoAcid.isSideChainAtom(donor);
    }

    @Override
    public double getEnergyContribution() {
        return 25;
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
