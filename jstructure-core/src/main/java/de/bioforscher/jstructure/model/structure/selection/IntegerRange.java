package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.model.Pair;

/**
 * Represents a defined range of residues identified by their residue numbers. Both indices are inclusive.
 * Created by S on 21.11.2016.
 */
public class IntegerRange extends Pair<Integer, Integer> {
    public IntegerRange(int i1, int i2) {
        super(i1, i2);
    }

    @Override
    public String toString() {
        return "Range [" + getLeft() + ", " + getRight() + "]";
    }
}