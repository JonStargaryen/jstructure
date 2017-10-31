package de.bioforscher.jstructure.mathematics.graph.partitioning;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a module or cluster or partition or community in a graph or network. Basically just a association of
 * amino acids and the contract, that for a chain each amino acid can only be assigned to one module. As modules are
 * basically subgraphs of a parent graph, this is reflected by this data structure which is derived from the standard
 * {@link Graph} implementation.
 */
public class Module<N> extends Graph<N> {
    private final Graph<N> parentGraph;
    private final String id;

    public Module(String id,
                  Graph<N> parentGraph,
                  List<N> nodes) {
        super(new Graph<>(nodes, parentGraph.getEdges()
                .stream()
                .filter(edge -> nodes.contains(edge.getLeft()) && nodes.contains(edge.getRight()))
                .collect(Collectors.toList())));
        this.parentGraph = parentGraph;
        this.id = id;
    }

    public String getIdentifier() {
        return id;
    }

    public Graph<N> getParentGraph() {
        return parentGraph;
    }

    public String getId() {
        return id;
    }

    /**
     * The number of edges to nodes of this module.
     * @param node the node to evaluate
     * @return the number of edges of a given node within this module
     */
    public int getIntraModuleDegree(N node) {
//        return (int) getEdges()
//                .stream()
//                .filter(edge -> edge.contains(node))
//                .count();
        return getDegreeOf(node);
    }

    /**
     * The number of edges to nodes of other modules.
     * @param node the node to evaluate
     * @return the number of edges of a given node to other modules
     */
    public int getInterModuleDegree(N node) {
        return (int) getInterModuleEdges()
                .stream()
                .filter(edge -> edge.contains(node))
                .count();
    }

    public List<Edge<N>> getIntraModuleEdges() {
        return getEdges();
    }

    public List<Edge<N>> getInterModuleEdges() {
        return parentGraph.getEdges()
                .stream()
                .filter(edge -> !getEdges().contains(edge))
                .collect(Collectors.toList());
    }
}
