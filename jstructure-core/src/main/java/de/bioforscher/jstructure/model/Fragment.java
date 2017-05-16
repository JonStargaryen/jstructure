package de.bioforscher.jstructure.model;

import java.util.List;
import java.util.stream.Collectors;

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

    @Override
    public String toString() {
        return elements.stream().map(Object::toString).collect(Collectors.joining(", ", "[", "]"));
    }
}