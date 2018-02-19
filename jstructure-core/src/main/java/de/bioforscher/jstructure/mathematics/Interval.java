package de.bioforscher.jstructure.mathematics;

public class Interval<N extends Number> extends Pair<N, N> {
    public Interval(N start, N endInclusive) {
        super(start, endInclusive);
    }

    @Override
    public String toString() {
        return "Interval [" + getLeft() + ", " + getRight() + "]";
    }
}
