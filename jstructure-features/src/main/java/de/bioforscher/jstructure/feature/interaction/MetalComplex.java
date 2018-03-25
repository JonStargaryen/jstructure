package de.bioforscher.jstructure.feature.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.nodes.Element;

import java.util.stream.Stream;

/**
 * A metal complex as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class MetalComplex extends PLIPInteraction {
    private Atom atom1;
    private Atom atom2;
    private double distance;
    private double rms;
    private String geometry;
    private String location;
    private String metalType;
    private int coordination;
    private int complexnum;

    MetalComplex(Group group, Element describingElement) {
        super(group, describingElement);
        this.atom1 = resolveAtom("metal_idx");
        this.metalType = getStringValueOfTag("metal_type");
        this.atom2 = resolveAtom("target_idx");
        this.coordination = getIntValueOfTag("coordination");
        this.distance = getDoubleValueOfTag("dist");
        this.location = getStringValueOfTag("location");
        this.rms = getDoubleValueOfTag("rms");
        this.geometry = getStringValueOfTag("geometry");
        this.complexnum = getIntValueOfTag("complexnum");
        this.partner2 = atom1.getParentGroup().equals(partner1) ? atom2.getParentGroup() : atom1.getParentGroup();
    }

    public Atom getAtom1() {
        return atom1;
    }

    public Atom getAtom2() {
        return atom2;
    }

    public double getDistance() {
        return distance;
    }

    public double getRms() {
        return rms;
    }

    public String getGeometry() {
        return geometry;
    }

    public String getLocation() {
        return location;
    }

    public String getMetalType() {
        return metalType;
    }

    public int getCoordination() {
        return coordination;
    }

    public int getComplexnum() {
        return complexnum;
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
