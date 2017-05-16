package de.bioforscher.jstructure.parser.plip.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

/**
 * A metal complex as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class MetalComplex extends PLIPInteraction {
    private Atom atom1, atom2;
    private double distance, rms;
    private String geometry, location, metalType;
    private int coordination, complexnum;

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
    public String toString() {
        return getClass().getSimpleName() + " atom1=" + toString(atom1) + " atom2=" + toString(atom2);
    }
}
