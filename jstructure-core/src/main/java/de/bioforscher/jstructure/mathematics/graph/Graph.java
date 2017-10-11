package de.bioforscher.jstructure.mathematics.graph;

import java.util.List;

/**
 * A really basic implementation of a graph.
 * @param <N> the type of nodes this graph contains
 */
public class Graph<N> {
    private final List<N> nodes;
    private final List<Edge<N>> edges;

    public Graph(Graph<N> graph) {
        this.nodes = graph.nodes;
        this.edges = graph.edges;
    }

    public Graph(List<N> nodes, List<Edge<N>> edges) {
        this.nodes = nodes;
        this.edges = edges;
    }

    public boolean containsNode(N node) {
        return nodes.contains(node);
    }

    public boolean containsEdge(Edge<N> edge) {
        return edges.contains(edge);
    }

    public int getNumberOfNodes() {
        return nodes.size();
    }

    public List<N> getNodes() {
        return nodes;
    }

    public List<Edge<N>> getEdges() {
        return edges;
    }
}
