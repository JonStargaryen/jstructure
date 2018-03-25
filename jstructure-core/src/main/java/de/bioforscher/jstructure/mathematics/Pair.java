package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.feature.FeatureContainer;

/**
 * Implementation of a unordered pair of objects. No ordering is assumed, however, as the entries can be of
 * different content types, one is referred to as left and the other as right entry.
 * @param <L> the content type of the left element
 * @param <R> the content type of the right element
 * Created by S on 02.10.2016.
 */
public class Pair<L, R> extends FeatureContainer {
    private final L left;
    private final R right;

    public Pair(L t1, R t2) {
        this.left = t1;
        this.right = t2;
    }

    public L getLeft() {
        return left;
    }

    public R getRight() {
        return right;
    }

    /**
     * Returns a new pair container with reversed order.
     * @return a new pair where the first entry is now second and vice versa
     */
    public Pair<R, L> flip() {
        Pair<R, L> flippedPair = new Pair<>(right, left);
        flippedPair.setFeatures(getFeatures());
        return flippedPair;
    }

    /**
     * Tests whether an object is an entry of this pair by employing {@link Object#equals(Object)}.
     * @param object the instance to check
     * @return true if either of this pair's entries equals the given object
     */
    public boolean contains(Object object) {
        return left.equals(object) || right.equals(object);
    }

    @Override
    public String toString() {
        return "(" + left + ", " + right + ")";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Pair<?, ?> pair = (Pair<?, ?>) o;

        return (left != null ? left.equals(pair.left) : pair.left == null) && (right != null ? right.equals(pair.right) : pair.right == null);
    }

    @Override
    public int hashCode() {
        int result = left != null ? left.hashCode() : 0;
        result = 31 * result + (right != null ? right.hashCode() : 0);
        return result;
    }
}