package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.interfaces.ShortestPathAlgorithm;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultEdge;

import java.util.*;
import java.util.stream.Collectors;

public class ProteinGraphCalculations {
    private final ProteinGraph graph;
    private final List<AminoAcid> nodes;
    private final Map<AminoAcid, ShortestPathAlgorithm.SingleSourcePaths<AminoAcid, DefaultEdge>> shortestPaths;
    private final DijkstraShortestPath<AminoAcid, DefaultEdge> dijkstraShortestPath;
    private final int numberOfNodes;
    private final double numberOfNodePairs;

    public ProteinGraphCalculations(ProteinGraph graph) {
        this.graph = graph;
        this.nodes = new ArrayList<>(graph.vertexSet());
        this.shortestPaths = new HashMap<>();
        this.dijkstraShortestPath = new DijkstraShortestPath<>(graph);
        for(AminoAcid aminoAcid : graph.vertexSet()) {
            this.shortestPaths.put(aminoAcid, dijkstraShortestPath.getPaths(aminoAcid));
        }
        this.numberOfNodes = graph.vertexSet().size();
        this.numberOfNodePairs = numberOfNodes * (numberOfNodes - 1) * 0.5;
    }

    /**
     * Determines the average path length between all pairs of nodes in a given graph. Graph must be connected.
     * See Vendruscolo, 2002 for definition.
     * @return the average path length of this graph
     */
    public double averageGraphPathLength() {
        double np = numberOfNodes * (numberOfNodes - 1) * 0.5;
        return 1 / np * graph.vertexSet()
                .stream()
                .mapToDouble(this::averageGraphPathLength)
                .average()
                .orElseThrow(() -> new IllegalArgumentException("could not compute average path length"));
    }

    public double averageGraphPathLength(AminoAcid source) {
        return graph.vertexSet()
                .stream()
                .mapToInt(target -> determineShortestPath(source, target).getLength())
                .average()
                .orElse(0);
    }

    /**
     * Betweenness is defined as the number shortest paths on the graph passing through this node normalized by the
     * total number of pairs of nodes of the graph.
     * See Vendruscolo, 2002 for definition.
     * @param node the node to evaluate
     * @return the betweenness of this node
     */
    public double betweenness(AminoAcid node) {
        return  shortestPathsPassingThrough(node).size() / numberOfNodePairs;
    }

    public double betweenness(ResidueIdentifier residueIdentifier) {
        return betweenness(resolve(residueIdentifier));
    }

    private AminoAcid resolve(ResidueIdentifier residueIdentifier) {
        return nodes.stream()
                .filter(node -> node.getResidueIdentifier().equals(residueIdentifier))
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("did not find residue with id " + residueIdentifier));
    }

    /**
     * The maximal path length from this node to any other node within the graph.
     * See Amitai, 2004 for definition.
     * @param node the node to evaluate
     * @return the maximal path length involving this node
     */
    public double closeness(AminoAcid node) {
        return 1 / graph.vertexSet().stream()
                .map(target -> determineShortestPath(node, target))
                .mapToInt(GraphPath::getLength)
                .average()
                .orElseThrow(() ->  new IllegalArgumentException("cannot evaluate closeness as graph is not fully connected"));
    }

    public double closeness(ResidueIdentifier residueIdentifier) {
        return closeness(resolve(residueIdentifier));
    }

    /**
     * Evaluates all nodes of a given graph by comparing the actual number of edges in the network of direct
     * neighbors to the theoretical maximum.
     * See Vendruscolo, 2002 for definition.
     * @return the clustering coefficient of this graph
     */
    public double clusteringCoefficient() {
        return graph.vertexSet().stream()
                .mapToDouble(this::clusteringCoefficient)
                .average()
                .orElseThrow(() -> new IllegalArgumentException("could not compute clustering coefficient of graph " +
                        graph + " as the graph does not contain edges"));
    }

    /**
     * The local clustering coefficient. For a given node and all its neighbors: how many edges are present compared
     * to the theoretical limit?
     * @param node the node to evaluate
     * @return the fraction of edges present
     */
    public double clusteringCoefficient(AminoAcid node) {
        List<AminoAcid> neighbors = getNeighborsFor(node);
        int numberOfNeighbors = neighbors.size();
        if(numberOfNeighbors == 1) {
            return 0;
        }
        long actualNumberOfEdges = SetOperations.uniquePairsOf(neighbors)
                .filter(pair -> containsEdge(pair.getLeft(), pair.getRight()))
                .count();
        return actualNumberOfEdges / (0.5 * (numberOfNeighbors * (numberOfNeighbors - 1)));
    }

    public double clusteringCoefficient(ResidueIdentifier residueIdentifier) {
        return clusteringCoefficient(resolve(residueIdentifier));
    }

    public int distinctNeighborhoodCount(AminoAcid aminoAcid) {
        List<Integer> interactingGroups = graph.getContacts().stream()
                .filter(plipInteraction -> aminoAcid.equals(plipInteraction.getLeft()) || aminoAcid.equals(plipInteraction.getRight()))
                .filter(inter -> Math.abs(inter.getLeft().getResidueIndex() - inter.getRight().getResidueIndex()) > 5)
                .map(plipInteraction -> aminoAcid.equals(plipInteraction.getLeft()) ? plipInteraction.getRight() : plipInteraction.getLeft())
                .map(Group::getResidueIndex)
                .collect(Collectors.toList());
        List<Integer> distinctNeighborhoods = new ArrayList<>();
        for(int interactingGroup : interactingGroups) {
            if(distinctNeighborhoods.stream()
                    .noneMatch(residueIndex -> Math.abs(residueIndex - interactingGroup) > 5)) {
                distinctNeighborhoods.add(interactingGroup);
            }
        }
        return distinctNeighborhoods.size();
    }

    public int distinctNeighborhoodCount(ResidueIdentifier residueIdentifier) {
        return distinctNeighborhoodCount(resolve(residueIdentifier));
    }

    private List<GraphPath<AminoAcid, DefaultEdge>> shortestPathsPassingThrough(AminoAcid node) {
        return SetOperations.uniquePairsOf(nodes)
                .map(pair -> shortestPaths.get(pair.getLeft()).getPath(pair.getRight()))
                .filter(graphPath -> graphPath.getVertexList().contains(node))
                .collect(Collectors.toList());
    }

    private List<AminoAcid> getNeighborsFor(AminoAcid node) {
        return graph.vertexSet()
                .stream()
                .filter(otherNode -> containsEdge(node, otherNode))
                .collect(Collectors.toList());
    }

    private boolean containsEdge(AminoAcid v1, AminoAcid v2) {
        return graph.containsEdge(v1, v2) || graph.containsEdge(v2, v1);
    }

    private GraphPath<AminoAcid, DefaultEdge> determineShortestPath(AminoAcid source, AminoAcid target) {
        try {
            return shortestPaths.get(source).getPath(target);
        } catch (NullPointerException e) {
            e.printStackTrace();
            return null;
        }
    }
}
