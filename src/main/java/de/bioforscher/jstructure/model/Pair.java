package de.bioforscher.jstructure.model;

import java.util.List;
import java.util.ArrayList;
import java.util.Collection;
import java.util.stream.Stream;

/**
 * Implementation of a unordered pair of equal objects.
 * Created by S on 02.10.2016.
 */
public class Pair<T> {
    private final T t1;
    private final T t2;

    public Pair(T t1, T t2) {
        this.t1 = t1;
        this.t2 = t2;
    }

    public T getFirst() {
        return t1;
    }

    public T getSecond() {
        return t2;
    }

    @Override
    public String toString() {
        return "(" + t1 + ", " + t2 + ")";
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
//    public static <T> Stream<Pair<T>> unorderPairsOf(List<T> collection) {
//        return IntStream.range(0, collection.size()).mapToObj(index1 -> IntStream.range(index1 + 1,
//                collection.size() + 1).mapToObj(index2 -> new Pair<>(collection.get(index1),
//                collection.get(index2))));
//    }
    public static <T> Stream<Pair<T>> unorderedPairsOf(List<T> collection) {
        List<Pair<T>> unorderedPairs = new ArrayList<>();
        for(int i = 0; i < collection.size() - 1; i++) {
            for(int j = i + 1; j <  collection.size(); j++) {
                unorderedPairs.add(new Pair(collection.get(i), collection.get(j)));
            }
        }
        return unorderedPairs.stream();
    }

    /**
     * Composes the cartesian product of two sets. They have the same content type.
     * @param collection1 a set of elements
     * @param collection2 another set of elements
     * @param <T> the content of the collections
     * @return the cartesian product of both sets
     */
    public static <T> Stream<Pair<T>> cartesianProductOf(Collection<T> collection1, Collection<T> collection2) {
        return collection1.stream().flatMap(element1 ->
            collection2.stream().map(element2 -> new Pair<>(element1, element2))
        );
    }
}