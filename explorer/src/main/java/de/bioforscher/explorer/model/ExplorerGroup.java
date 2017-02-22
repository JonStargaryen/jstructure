package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The reduced representation of {@link Group} objects.
 * Created by bittrich on 2/22/17.
 */
public class ExplorerGroup {
    private int resn;
    private List<ExplorerAtom> atoms;
    private String olc, tlc, type, sse;
    private double rasa;

    public ExplorerGroup() {

    }

    public ExplorerGroup(Group group) {
        this.resn = group.getResidueNumber();
        this.atoms = group.atoms()
                .map(ExplorerAtom::new)
                .collect(Collectors.toList());
        this.olc = group.getGroupInformation().getOneLetterCode();
        this.tlc = group.getGroupInformation().getThreeLetterCode();

        if(group.isAminoAcid()) {
            this.type = "aa";
            // assign features
            this.rasa = group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA);
            this.sse = group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().getOneLetterRepresentation();
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

    public String getSse() {
        return sse;
    }

    public double getRasa() {
        return rasa;
    }
}
