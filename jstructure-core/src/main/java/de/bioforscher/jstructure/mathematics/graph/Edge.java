package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.model.Pair;

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
}
