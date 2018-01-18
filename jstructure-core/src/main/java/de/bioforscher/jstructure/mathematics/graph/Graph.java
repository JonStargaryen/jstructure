package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.mathematics.Calculable;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A really basic implementation of a graph.
 * @param <N> the type of nodes this graph contains
 */
public class Graph<N> implements Calculable<LinearAlgebra.GraphLinearAlgebra> {
    private final List<N> nodes;
    private final List<Edge<N>> edges;

    public Graph(Graph<N> graph) {
        this.nodes = graph.nodes;
        this.edges = graph.edges;
    }

    public Graph(List<N> nodes) {
        this.nodes = nodes;
        this.edges = new ArrayList<>();
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

    public boolean containsEdge(Pair<N, N> pair) {
        return containsEdge(pair.getLeft(), pair.getRight());
    }

    public boolean containsEdge(N node1, N node2) {
        return edges()
                .anyMatch(edge -> edge.contains(node1) && edge.contains(node2));
    }

    public int getNumberOfNodes() {
        return nodes.size();
    }

    public int getNumberOfEdges() {
        return edges.size();
    }

    public Stream<N> nodes() {
        return nodes.stream();
    }

    public List<N> getNodes() {
        return nodes;
    }

    public Stream<Edge<N>> edges() {
        return edges.stream();
    }

    public List<Edge<N>> getEdges() {
        return edges;
    }

    public Stream<Edge<N>> edgesFor(N node) {
        return edges()
                .filter(edge -> edge.contains(node));
    }

    public int getDegreeOf(N node) {
        return (int) edgesFor(node)
                .count();
    }

    public List<Edge<N>> getEdgesFor(N node) {
        return edgesFor(node)
                .collect(Collectors.toList());
    }

    public Stream<N> neighborsFor(N node) {
        return nodes()
                .filter(candiate -> !node.equals(candiate))
                .filter(candidate -> edgeBetween(candidate, node).isPresent());
    }

    public List<N> getNeighborsFor(N node) {
        return neighborsFor(node)
                .collect(Collectors.toList());
    }

    public Stream<Pair<Edge<N>, Edge<N>>> commonNeighborhood(Edge<N> edge) {
        N node1 = edge.getLeft();
        N node2 = edge.getRight();

        return nodes()
                .filter(node -> !node.equals(node1) && !node.equals(node2))
                .map(node -> new Pair<>(edgeBetween(node, node1), edgeBetween(node, node2)))
                .filter(pair -> pair.getLeft().isPresent() && pair.getRight().isPresent())
                .map(pair -> new Pair<>(pair.getLeft().get(), pair.getRight().get()));
    }

    public int getSizeOfCommonNeighborhood(Edge<N> edge) {
        return (int) commonNeighborhood(edge)
                .count();
    }

    private Optional<Edge<N>> edgeBetween(N n1, N n2) {
        return edges()
                .filter(edge -> edge.contains(n1) && edge.contains(n2))
                .findFirst();
    }


    public void remove(N nodeToRemove) {
        nodes.remove(nodeToRemove);
    }

    public void remove(Edge<N> edgeToRemove) {
        edges.remove(edgeToRemove);
    }

    @Override
    public LinearAlgebra.GraphLinearAlgebra<N> calculate() {
        return LinearAlgebra.on(this);
    }
}
