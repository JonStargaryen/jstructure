package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.model.structure.Group;

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
        this.resn = group.getResidueNumber();
        this.aa = group.getGroupInformation().getOneLetterCode();
        this.rasa = DECIMAL_FORMAT.format(group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA));
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
