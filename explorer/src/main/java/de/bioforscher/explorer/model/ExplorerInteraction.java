package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

/**
 * PLIP data.
 * Created by bittrich on 4/6/17.
 */
public class ExplorerInteraction {
    private int partner1, partner2;
    private double[] coords1, coords2, center;
    private boolean isLongRange;

    public ExplorerInteraction() {
    }

    public ExplorerInteraction(PLIPInteraction interaction, LinearAlgebraAtom.Transformation transformation) {
        this.partner1 = interaction.getPartner1().getResidueNumber();
        this.partner2 = interaction.getPartner2().getResidueNumber();
        //TODO remove
        if (transformation != null) {
            this.coords1 = LinearAlgebra3D.add(LinearAlgebra3D.multiply(interaction.getCoords1(), transformation.getRotation()), transformation.getTranslation());
            this.coords2 = LinearAlgebra3D.add(LinearAlgebra3D.multiply(interaction.getCoords2(), transformation.getRotation()), transformation.getTranslation());
            this.center = LinearAlgebra3D.add(LinearAlgebra3D.multiply(interaction.getRepresentation(), transformation.getRotation()), transformation.getTranslation());
        } else {
            this.coords1 = interaction.getCoords1();
            this.coords2 = interaction.getCoords2();
            this.center = interaction.getRepresentation();
        }
        this.isLongRange = Math.abs(partner1 - partner2) > 5;
    }

    public int getPartner1() {
        return partner1;
    }

    public int getPartner2() {
        return partner2;
    }

    public double[] getCoords1() {
        return coords1;
    }

    public double[] getCoords2() {
        return coords2;
    }

    public double[] getCenter() {
        return center;
    }

    public boolean isLongRange() {
        return isLongRange;
    }
}
