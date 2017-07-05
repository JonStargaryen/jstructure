package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.nodes.Element;

import java.util.stream.Stream;

/**
 * A water bridge as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class WaterBridge extends PLIPInteraction {
    private Atom atom1;
    private Atom atom2;
    private double distanceAW;
    private double distanceDW;
    private double donorAngle;
    private double waterAngle;
    private boolean protIsDon;

    WaterBridge(Group group, Element describingElement) {
        super(group, describingElement);
        this.distanceAW = getDoubleValueOfTag("dist_a-w");
        this.distanceDW = getDoubleValueOfTag("dist_d-w");
        this.donorAngle = getDoubleValueOfTag("don_angle");
        this.waterAngle = getDoubleValueOfTag("water_angle");
        this.protIsDon = getBooleanValueOfTag("protisdon");
        this.atom1 = resolveAtom("donor_idx");
        this.atom2 = resolveAtom("acceptor_idx");
        this.partner2 = protIsDon ? atom2.getParentGroup() : atom1.getParentGroup();
    }

    public Atom getAtom1() {
        return atom1;
    }

    public Atom getAtom2() {
        return atom2;
    }

    public double getDistanceAW() {
        return distanceAW;
    }

    public double getDistanceDW() {
        return distanceDW;
    }

    public double getDonorAngle() {
        return donorAngle;
    }

    public double getWaterAngle() {
        return waterAngle;
    }

    public boolean isProtIsDon() {
        return protIsDon;
    }

    @Override
    public boolean isBackboneInteraction() {
        return AminoAcid.isBackboneAtom(atom1) && AminoAcid.isBackboneAtom(atom2);
    }

    @Override
    public boolean isSideChainInteraction() {
        return AminoAcid.isSideChainAtom(atom1) && AminoAcid.isSideChainAtom(atom2);
    }

    @Override
    boolean isSane() {
        return atom1 != null && atom2 != null;
    }

    @Override
    public Stream<Atom> allAtoms() {
        return Stream.of(atom1, atom2);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " atom1=" + toString(atom1) + " atom2=" + toString(atom2);
    }
}
