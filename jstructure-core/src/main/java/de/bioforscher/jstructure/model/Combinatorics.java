package de.bioforscher.jstructure.model;

import java.util.Collection;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Access to combinatoric capabilities such as constructing {@link Fragment} or {@link Pair} objects.
 * Created by S on 15.11.2016.
 */
public class Combinatorics {
    private Combinatorics() {
        // deny instantiation
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
    public static <T> Stream<Pair<T, T>> uniquePairsOf(final List<T> collection) {
        return IntStream.range(0, collection.size() - 1)
                .boxed()
                .flatMap(i -> IntStream.range(i + 1, collection.size())
                        .mapToObj(j -> new Pair<>(collection.get(i), collection.get(j))));
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
        return collection1.stream()
                .flatMap(element1 -> collection2.stream()
                        .map(element2 -> new Pair<>(element1, element2)));
    }

    /**
     * Combines two collections into one select where each pair describes the sequential combination of each element at
     * the corresponding position of the input collections.
     *
     * <pre>
     *     A = { a1, a2, a3 }
     *     B = { b1, b2, b3 }
     *     A|B = { a1b1, a2b2, a3b3 }
     * </pre>
     *
     * Both collections must match in size.
     * @param collection1 a collection of elements
     * @param collection2 a similiarly size collections of other elements
     * @param <L> the type of the first collection
     * @param <R> the type of the second collection
     * @return all sequential pairs of the given collections
     */
    public static <L, R> Stream<Pair<L, R>> sequentialPairsOf(List<L> collection1, List<R> collection2) {
        if(collection1.size() != collection2.size()) {
            throw new IllegalArgumentException("collections must agree in size. found " + collection1.size() + " and " +
                    collection2.size());
        }

        return IntStream.range(0, collection1.size())
                .mapToObj(i -> new Pair<>(collection1.get(i), collection2.get(i)));
    }

    /**
     * Fragments a list into a select of sublists of a desired size.
     * <pre>
     *     A = { a, b, c, d }
     *     will be fragmented into
     *     {{ a, b }, { b, c }, { c, d }}
     *     for a fragmentSize of 2
     *     and into
     *     {{ a, b, c }, { b, c, d }}
     *     for a fragmentSize of 3.
     * </pre>
     * @param elements a collection of elements
     * @param fragmentSize the desired number of grouped elements
     * @param <T> the content type
     * @return a select of fragments which can be composed based on the given collection
     */
    public static <T> Stream<Fragment<T>> fragmentsOf(List<T> elements, int fragmentSize) {
        if(fragmentSize < 2) {
            throw new IllegalArgumentException("fragment size cannot be smaller than 2 - found " + fragmentSize);
        }
        if(fragmentSize >= elements.size()) {
            throw new IllegalArgumentException("fragment size cannot be equal or exceed number of elements - found " +
                    fragmentSize + " and " + elements.size() + " elements");
        }

        return IntStream.range(0, elements.size() - fragmentSize + 1)
                .mapToObj(i -> new Fragment<>(elements.subList(i, i + fragmentSize)));
    }
}
