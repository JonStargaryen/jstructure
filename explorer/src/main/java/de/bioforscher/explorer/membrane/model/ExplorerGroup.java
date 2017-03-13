package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.sse.DSSPSecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Group} objects.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerGroup {
    private int resn;
    private List<ExplorerAtom> atoms;
    private String olc, tlc, type;
    private double rasa, sse;

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
                this.rasa = group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA);
            } catch (NullPointerException e) {
                this.rasa = 0.0;
            }
            try {
                this.sse = group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().ordinal() / (double) DSSPSecondaryStructureElement.values().length;
            } catch (NullPointerException e) {
                this.sse = 0.0;
            }
        } else if(group.isNucleotide()) {
            this.type = "nucl";
        } else {
            this.type = "het";
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

    public double getSse() {
        return sse;
    }

    public double getRasa() {
        return rasa;
    }
}
