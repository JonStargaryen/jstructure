package de.bioforscher.jstructure.model;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Implementation of an ordered consecutive list of equal objects.
 * Created by S on 02.10.2016.
 */
public class Fragment<T> {
    private List<T> elements;

    public Fragment(List<T> elements) {
        this.elements = elements;
    }

    /**
     *
     * @param index
     * @return
     */
    public T getElement(int index) {
        return elements.get(index);
    }

    @Override
    public String toString() {
        return elements.stream().map(Object::toString).collect(Collectors.joining(", ", "[", "]"));
    }

    /**
     * Fragments a list into a stream of sublists of a desired size.
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
     * @return a stream of fragments which can be composed based on the given collection
     */
    public static <T> Stream<Fragment<T>> fragmentsOf(List<T> elements, int fragmentSize) {
        if(fragmentSize < 2) {
            throw new IllegalArgumentException("fragment size cannot be smaller than 2 - found " + fragmentSize);
        }
        if(fragmentSize >= elements.size()) {
            throw new IllegalArgumentException("fragment size cannot be equal or exceed number of elements - found " +
                    fragmentSize + " and " + elements.size() + " elements");
        }

        return IntStream.range(0, elements.size() - fragmentSize + 1).mapToObj(i -> new Fragment<>(elements.subList(i,
                i + fragmentSize)));
    }
}