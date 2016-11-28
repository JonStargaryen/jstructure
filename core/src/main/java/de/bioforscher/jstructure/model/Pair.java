package de.bioforscher.jstructure.model;

/**
 * Implementation of a unordered pair of equal objects. No ordering is assumed, however, as the entries can be of
 * different content types, one is referred to as left and the other as right entry.
 * @param <L> the content type of the left element
 * @param <R> the content type of the right element
 * Created by S on 02.10.2016.
 */
public class Pair<L, R> {
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
        return new Pair<>(right, left);
    }

    /**
     * Tests whether an object is an entry of this pair.
     * @param object the instance to check
     * @return true if either of this pair's entries equals the given object
     */
    public boolean contains(Object object) {
        //TODO contract which Object#equals() is employed, does the reference have to be identical or mere content?
        return left.equals(object) || right.equals(object);
    }

    @Override
    public String toString() {
        return "(" + left + ", " + right + ")";
    }
}