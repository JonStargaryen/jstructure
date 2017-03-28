package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

/**
 * PLIP data.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerInteraction {
    private double[] coords1, coords2;

    public ExplorerInteraction() {
    }

    public ExplorerInteraction(PLIPInteraction interaction, LinearAlgebraAtom.Transformation transformation) {
        this.coords1 = LinearAlgebra3D.add(LinearAlgebra3D.multiply(interaction.getCoords1(), transformation.getRotation()), transformation.getTranslation());
        this.coords2 = LinearAlgebra3D.add(LinearAlgebra3D.multiply(interaction.getCoords2(), transformation.getRotation()), transformation.getTranslation());
    }

    public double[] getCoords1() {
        return coords1;
    }

    public double[] getCoords2() {
        return coords2;
    }
}
