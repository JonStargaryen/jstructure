package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

/**
 * PLIP data.
 * Created by bittrich on 3/20/17.
 */
public class ExplorerInteraction {
    private double[] coords1, coords2;

    public ExplorerInteraction() {
    }

    public ExplorerInteraction(PLIPInteraction interaction) {
        this.coords1 = interaction.getCoords1();
        this.coords2 = interaction.getCoords2();
    }

    public double[] getCoords1() {
        return coords1;
    }

    public double[] getCoords2() {
        return coords2;
    }
}
