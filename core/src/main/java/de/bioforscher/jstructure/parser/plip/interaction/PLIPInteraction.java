package de.bioforscher.jstructure.parser.plip.interaction;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.jsoup.nodes.Element;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * The abstract implementation of a PLIP interaction.
 * Created by bittrich on 2/15/17.
 */
public abstract class PLIPInteraction {
    private static final String TRUE_FLAG = "True";
    private Protein protein;
    Group partner1, partner2;
    private Element describingElement;
    private double[] coords1, coords2;

    /**
     * Enforce the signature of the constructor of all subclasses.
     *
     * @param group             one partner that is described directly by this element - could be parsed from the 2nd argument, then
     *                          however we would be stuck in a constructor of an object which cannot be created
     * @param describingElement the xml element to parse to gather all information
     */
    PLIPInteraction(Group group, Element describingElement) {
        this.partner1 = group;
        this.describingElement = describingElement;
        this.protein = partner1.getParentChain().getParentProtein();
        Element ligcoo = describingElement.getElementsByTag("ligcoo").first();
        this.coords1 = extractCoordinateArray(ligcoo);
        Element protcoo = describingElement.getElementsByTag("protcoo").first();
        this.coords2 = extractCoordinateArray(protcoo);
    }

    private double[] extractCoordinateArray(Element coo) {
        return new double[] { Double.valueOf(coo.child(0).text()), Double.valueOf(coo.child(1).text()), Double.valueOf(coo.child(2).text()) };
    }


    String getStringValueOfTag(String tagName) {
        return describingElement.getElementsByTag(tagName).first().text();
    }

    double getDoubleValueOfTag(String tagName) {
        return Double.valueOf(getStringValueOfTag(tagName));
    }

    int getIntValueOfTag(String tagName) {
        return Integer.valueOf(getStringValueOfTag(tagName));
    }

    boolean getBooleanValueOfTag(String tagName) {
        return getStringValueOfTag(tagName).equals(TRUE_FLAG);
    }

    Atom resolveAtom(String tagName) {
        int pdbSerial = getIntValueOfTag(tagName);
        Optional<Atom> atom = Selection.on(protein)
                .pdbSerial(pdbSerial)
                .asOptionalAtom();

        if(atom.isPresent()) {
            return atom.get();
        }

        throw new NoSuchElementException("no atom found for pdbSerial '" + pdbSerial + "'");
    }

    List<Atom> resolveAtoms(String tagname) {
        return describingElement.getElementsByTag(tagname).stream()
                .map(Element::text)
                .map(Integer::valueOf)
                .map(atomNumber -> Selection.on(protein)
                        .pdbSerial(atomNumber)
                        .asOptionalAtom())
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }

    public Group getPartner1() {
        return partner1;
    }

    public Group getPartner2() {
        return partner2;
    }

    /**
     * True iff the interaction is observed within the same group.
     * @return true when <code>partner1.equals(partner2)</code>
     */
    public boolean partnerIsInSameGroup() {
        return partner1.equals(partner2);
    }

    /**
     * True iff this is an interaction to another amino acid.
     * @return true when <code>partner2.isAminoAcid()</code> evaluates to <code>true</code>
     */
    public boolean partnerIsAminoAcid() {
        return partner2.isAminoAcid();
    }

    String toString(Atom atom) {
        return "'" + atom.getParentGroup().getIdentifier() + "-" + atom.getIdentifier() + "'";
    }

    String toString(List<Atom> atoms) {
        return "'" + atoms.stream()
                .map(atom -> atom.getParentGroup().getIdentifier() + "-" + atom.getIdentifier())
                .collect(Collectors.joining(", ")) + "'";
    }

    public double[] getCoords1() {
        return coords1;
    }

    public double[] getCoords2() {
        return coords2;
    }
}
