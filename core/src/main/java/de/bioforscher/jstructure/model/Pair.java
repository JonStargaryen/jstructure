package de.bioforscher.jstructure.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

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

    @Override
    public String toString() {
        return "(" + left + ", " + right + ")";
    }

    /**
     * Creates a collection of all unordered pairs which can be created using the method's argument. An element will
     * never be paired with itself.
     * <pre>
     *     A = { a, b, c }
     *     [A]² = { ab, ac, bc }
     *     'ab' is equal to 'ba'
     *     'aa' is no valid combination produced by this method
     * </pre>
     * @param collection the input collection
     * @param <T> the content of the collection
     * @return all {@link Pair} objects which can be created from the collection
     */
    public static <T> Stream<Pair<T, T>> uniquePairsOf(List<T> collection) {
        //TODO streamify
        //TODO rename
//        IntStream.range(0, collection.size() - 1)
//                 .flatMap(i -> IntStream.range(i + 1, collection.size())
//                                        .mapToObj());

        List<Pair<T, T>> unorderedPairs = new ArrayList<>();
        for(int i = 0; i < collection.size() - 1; i++) {
            for(int j = i + 1; j <  collection.size(); j++) {
                unorderedPairs.add(new Pair<>(collection.get(i), collection.get(j)));
            }
        }
        return unorderedPairs.stream();
    }

    /**
     * Composes the cartesian product of two sets
     * @param collection1 a set of elements
     * @param collection2 another set of elements
     * @param <L> the content of the first collection
     * @param <R> the content of the second collection
     * @return the cartesian product of both sets
     */
    public static <L, R> Stream<Pair<L, R>> cartesianProductOf(Collection<L> collection1, Collection<R> collection2) {
        return collection1.stream().flatMap(element1 ->
            collection2.stream().map(element2 -> new Pair<>(element1, element2))
        );
    }

    public static <L, R> Stream<Pair<L, R>> sequentialPairsOf(List<L> collection1, List<R> collection2) {
        if(collection1.size() != collection2.size()) {
            throw new IllegalArgumentException("collections must agree in size. found " + collection1.size() + " and " +
                collection2.size());
        }

        return IntStream.range(0, collection1.size()).mapToObj(i -> new Pair<>(collection1.get(i), collection2.get(i)));
    }
}