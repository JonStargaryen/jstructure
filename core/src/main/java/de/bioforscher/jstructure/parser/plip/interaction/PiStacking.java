package de.bioforscher.jstructure.parser.plip.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

import java.util.ArrayList;
import java.util.List;

/**
 * A pi-stacking interaction as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class PiStacking extends PLIPInteraction {
    private List<Atom> atoms1, atoms2;
    private double distance, angle, offset;
    private String type;

    PiStacking(Group group, Element describingElement) {
        super(group, describingElement);
        this.distance = getDoubleValueOfTag("centdist");
        this.angle = getDoubleValueOfTag("angle");
        this.offset = getDoubleValueOfTag("offset");
        this.type = getStringValueOfTag("type");
        this.atoms1 = new ArrayList<>();
        this.atoms2 = resolveAtoms("idx");
        this.partner2 = atoms2.get(0).getParentGroup();
    }

    public List<Atom> getAtoms1() {
        return atoms1;
    }

    public List<Atom> getAtoms2() {
        return atoms2;
    }

    public double getDistance() {
        return distance;
    }

    public double getAngle() {
        return angle;
    }

    public double getOffset() {
        return offset;
    }

    public String getType() {
        return type;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " atoms1=" + toString(atoms1) + " atoms2=" + toString(atoms2);
    }
}
