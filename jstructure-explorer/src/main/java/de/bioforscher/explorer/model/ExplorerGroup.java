package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.text.DecimalFormat;

/**
 * Represents one group and all related information.
 * Created by bittrich on 4/20/17.
 */
public class ExplorerGroup {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####");

    private int resn;
    private String aa, rasa;

    public ExplorerGroup() {
    }

    public ExplorerGroup(Group group) {
        this.resn = group.getResidueIdentifier().getResidueNumber();
        this.aa = group.isAminoAcid() ? ((AminoAcid) group).getOneLetterCode() : group.getThreeLetterCode();
        try {
            this.rasa = DECIMAL_FORMAT.format(group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea());
        } catch (NullPointerException e) {
            // happens for mutation sites
        }
    }

    public int getResn() {
        return resn;
    }

    public String getAa() {
        return aa;
    }

    public String getRasa() {
        return rasa;
    }
}
