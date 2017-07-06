package de.bioforscher.jstructure.model;

import java.util.List;
import java.util.stream.Stream;

/**
 * Implementation of an ordered consecutive collection of equal objects. Delegates to an internal list.
 * Created by S on 02.10.2016.
 */
public class Fragment<T> {
    private List<T> elements;

    public Fragment(List<T> elements) {
        this.elements = elements;
    }

    /**
     * Returns a specific element.
     * @param index the requested index
     * @return the specified element
     */
    public T getElement(int index) {
        return elements.get(index);
    }

    public List<T> getElements() {
        return elements;
    }

    public Stream<T> elements() {
        return elements.stream();
    }

    public int size() {
        return elements.size();
    }

    @Override
    public String toString() {
        return elements.toString();
    }
}