package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.mathematics.Pair;

import java.util.Objects;

public class Edge<N> extends Pair<N, N> {
    public static final double DEFAULT_WEIGHT = 1;

    private final double weight;

    public Edge(N n1, N n2) {
        this(n1, n2, DEFAULT_WEIGHT);
    }

    public Edge(N n1, N n2, double weight) {
        super(n1, n2);
        this.weight = weight;
    }

    public double getWeight() {
        return weight;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        Edge<?> edge = (Edge<?>) o;
        return Double.compare(edge.getWeight(), getWeight()) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), getWeight());
    }
}
