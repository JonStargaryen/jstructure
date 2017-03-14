package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.cerosene.SequenceCerosene;
import de.bioforscher.jstructure.feature.sse.DSSPSecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;

import static de.bioforscher.explorer.membrane.ColorMapping.*;

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
    private double[] rasa, sse, cerosene;

    public ExplorerGroup() {

    }

    ExplorerGroup(Group group) {
        this.resn = group.getResidueNumber();
        this.atoms = group.atoms()
                .map(ExplorerAtom::new)
                .collect(Collectors.toList());
        this.olc = group.getGroupInformation().getOneLetterCode();
        this.tlc = group.getGroupInformation().getThreeLetterCode();

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
}
