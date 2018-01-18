package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.mathematics.Fragment;

import java.util.List;

public class GraphPath<N> extends Fragment<N> {
    public GraphPath(List<N> elements) {
        super(elements);
    }

    /**
     * Reports the length of this path. Either the number of elements in the underlying collection or
     * {@link Integer#MAX_VALUE} if no path is defined.
     * @return the legnth of the described path
     */
    @Override
    public int size() {
        return getElements().isEmpty() ? Integer.MAX_VALUE : super.size();
    }
}
