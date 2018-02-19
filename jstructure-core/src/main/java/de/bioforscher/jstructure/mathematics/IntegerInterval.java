package de.bioforscher.jstructure.mathematics;

/**
 * Represents a defined range of residues identified by their residue numbers. Both indices are inclusive.
 * Created by S on 21.11.2016.
 */
public class IntegerInterval extends Interval<Integer> {
    public IntegerInterval(int i1, int i2) {
        super(i1, i2);
    }
}