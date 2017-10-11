package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.stream.Collectors;

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

    /**
     * Returns the NetCarto representation of this graph. I.e. a list of edges in the format:
     * <pre>
     *     nodeName1 nodeName2
     * </pre>
     * @return
     */
    public String toNetCartoString() {
        return edges.stream()
                .map(edge -> determineNodeName(edge.getLeft()) + " " + determineNodeName(edge.getRight()))
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private String determineNodeName(N node) {
        if(node instanceof AminoAcid) {
            return ((AminoAcid) node).getResidueIdentifier().toString();
        }
        return node.toString();
    }
}
