package de.bioforscher.jstructure.mathematics.graph.clustering;

import de.bioforscher.jstructure.mathematics.graph.Graph;

import java.util.ArrayList;

/**
 * Represents a module or cluster or partition or community in a graph or network. Basically just a association of
 * amino acids and the contract, that for a chain each amino acid can only be assigned to one module.
 */
public class Module<N> extends Graph<N> {
    private final String id;

    public Module(String id, Graph<N> subgraph) {
        super(subgraph);
        this.id = id;
        new ArrayList<>();
    }

    public String getIdentifier() {
        return id;
    }
}
