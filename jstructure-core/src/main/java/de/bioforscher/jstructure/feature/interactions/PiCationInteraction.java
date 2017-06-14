package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.nodes.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

/**
 * A pi-cation interaction as described by PLIP.
 * Created by bittrich on 2/15/17.
 */
public class PiCationInteraction extends PLIPInteraction {
    private List<Atom> atoms1;
    private List<Atom> atoms2;
    private double distance;
    private double offset;
    private boolean protcharged;
    private String ligandGroup;

    PiCationInteraction(Group group, Element describingElement) {
        super(group, describingElement);
        this.distance = getDoubleValueOfTag("dist");
        this.offset = getDoubleValueOfTag("offset");
        this.protcharged = getBooleanValueOfTag("protcharged");
        this.ligandGroup = getStringValueOfTag("lig_group");
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

    public double getOffset() {
        return offset;
    }

    public boolean isProtcharged() {
        return protcharged;
    }

    public String getLigandGroup() {
        return ligandGroup;
    }

    @Override
    public boolean isBackboneInteraction() {
        return Stream.concat(atoms1.stream(), atoms2.stream())
                .allMatch(AminoAcid::isBackboneAtom);
    }

    @Override
    public boolean isSideChainInteraction() {
        return Stream.concat(atoms1.stream(), atoms2.stream())
                .allMatch(AminoAcid::isSideChainAtom);
    }

    @Override
    boolean isSane() {
        return isSane(atoms1) && isSane(atoms2);
    }

    @Override
    Stream<Atom> allAtoms() {
        return Stream.concat(atoms1.stream(), atoms2.stream());
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " atoms1=" + toString(atoms1) + " atom2=" + toString(atoms2);
    }
}
