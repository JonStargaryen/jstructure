package de.bioforscher.explorer.membrane.o.dmodel.representative;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.cerosene.SequenceCerosene;
import de.bioforscher.jstructure.feature.sse.DSSPSecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;

import static de.bioforscher.explorer.membrane.o.dmodel.ColorMapping.*;

/**
 * The reduced representation of {@link Group} objects.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerGroup {
    private int resn;
    private List<ExplorerAtom> atoms;
    private String olc, tlc, type;
    /*
     * features to visualize are represented by double[] (i.e. hsv values)
     */
    private double[] rasa, sse, cerosene, aa;

    public ExplorerGroup() {

    }

    ExplorerGroup(Group group) {
        this.resn = group.getResidueNumber();
        this.atoms = group.atoms()
                .map(ExplorerAtom::new)
                .collect(Collectors.toList());
        this.olc = group.getGroupInformation().getOneLetterCode();
        this.tlc = group.getGroupInformation().getThreeLetterCode();
//        this.aa = discreteToHsv(AminoAcidFamily.valueOfIgnoreCase(olc).orElse(AminoAcidFamily.UNKNOWN).ordinal(), AminoAcidFamily.values().length);
        this.aa = mapAminoAcidToHsv(olc);

        if(group.isAminoAcid()) {
            this.type = "aa";
            // assign features - they are supposed to be normalized to the interval [0,1]
            try {
                this.rasa = continuousToHsv(group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA));
            } catch (NullPointerException e) {
                this.rasa = HSV_MISSING_VALUE;
            }
            try {
                this.sse = discreteToHsv(group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().ordinal(), DSSPSecondaryStructureElement.values().length);
            } catch (NullPointerException e) {
                this.sse = HSV_MISSING_VALUE;
            }
        } else if(group.isNucleotide()) {
            this.type = "nucl";
        } else {
            this.type = "het";
        }

        if(!this.type.equals("het")) {
            double[] hsv = group.getFeature(double[].class, SequenceCerosene.SEQUENCE_CEROSENE_REPRESENTATION);
            this.cerosene = new double[]{ hsv[0] * 360, hsv[1] * 100, hsv[2] * 100 };
        }
    }

    private double[] mapAminoAcidToHsv(String olc) {
        switch(olc) {
            case "H":
                return new double[] { 240, 100, 67 };
            case "R":
                return new double[] { 240, 100, 80 };
            case "N":case "Q":
                return new double[] { 105, 57, 47 };
            case "D":case "E":
                return new double[] { 98, 80, 67 };
            case "M":case "C":
                return new double[] { 60, 80, 50 };
            case "K":
                return new double[] { 207, 100, 70 };
            case "S":case "T":case "Y":
                return new double[] { 39, 100, 50 };
            case "A":case "I":case "L":case "F":case "P":case"V":case "W":case "G":
                return new double[] { 0, 0, 75 };
            default:
                return new double[]{ 206, 14, 22 };
        }
    }

    public int getResn() {
        return resn;
    }

    public List<ExplorerAtom> getAtoms() {
        return atoms;
    }

    public String getOlc() {
        return olc;
    }

    public String getTlc() {
        return tlc;
    }

    public String getType() {
        return type;
    }

    public double[] getSse() {
        return sse;
    }

    public double[] getRasa() {
        return rasa;
    }

    public double[] getCerosene() {
        return cerosene;
    }

    public double[] getAa() {
        return aa;
    }
}
