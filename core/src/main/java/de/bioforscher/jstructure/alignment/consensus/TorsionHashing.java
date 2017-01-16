package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.List;

/**
 * Created by bittrich on 1/9/17.
 */
public class TorsionHashing {
    public static final double DEFAULT_TORSION_ANGLE_BIN_SIZE = 10;
    private double torsionAngleBinSize;

    public TorsionHashing(double torsionAngleBinSize) {
        this.torsionAngleBinSize = torsionAngleBinSize;
    }

    public void composeClusterRepresentation(List<? extends AtomContainer> containers) {

    }
}
