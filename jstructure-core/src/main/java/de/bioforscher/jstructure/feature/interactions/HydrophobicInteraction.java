package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Element;

/**
 * A hydrophobic interaction as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class HydrophobicInteraction extends PLIPInteraction {
    private Atom atom1;
    private Atom atom2;
    private double distance;

    HydrophobicInteraction(Group group, Element describingElement) {
        super(group, describingElement);
        this.distance = getDoubleValueOfTag("dist");
        this.atom1 = resolveAtom("ligcarbonidx");
        this.atom2 = resolveAtom("protcarbonidx");
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

    @Override
    public String toString() {
        return getClass().getSimpleName() + " atom1=" + toString(atom1) + " atom2=" + toString(atom2);
    }
}
