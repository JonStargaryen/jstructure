package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.model.Pair;

/**
 * Represents a defined range of residues identified by their residue numbers. Both indices are inclusive.
 * Created by S on 21.11.2016.
 */
public class Range extends Pair<Integer, Integer> {
    public Range(int i1, int i2) {
        super(i1, i2);
    }
}