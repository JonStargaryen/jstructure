package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.interfaces.ShortestPathAlgorithm;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultEdge;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class TopologicPropertyCalculator extends FeatureProvider {
    @Override
    protected void processInternally(Structure structure) {
        structure.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        ProteinGraph fullPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015);
        ProteinGraph hydrogenPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROGEN_BONDS);
        ProteinGraph hydrophobicPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROPHOBIC_INTERACTION);
        ProteinGraph conventionalGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.CALPHA8);

        ProteinGraphCalculations fullPlipCalculations = new ProteinGraphCalculations(fullPlipGraph);
        ProteinGraphCalculations hydrogenPlipCalculations = new ProteinGraphCalculations(hydrogenPlipGraph);
        ProteinGraphCalculations hydrophobicPlipCalculations = new ProteinGraphCalculations(hydrophobicPlipGraph);
        ProteinGraphCalculations conventionalCalculations = new ProteinGraphCalculations(conventionalGraph);

        ChainTopologicProperties fullPlipChainTopologicProperties = new ChainTopologicProperties(this,
                fullPlipCalculations.averageGraphPathLength(),
                fullPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties hydrogenPlipChainTopologicProperties = new ChainTopologicProperties(this,
                hydrogenPlipCalculations.averageGraphPathLength(),
                hydrogenPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties hydrophobicPlipChainTopologicProperties = new ChainTopologicProperties(this,
                hydrophobicPlipCalculations.averageGraphPathLength(),
                hydrophobicPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties conventionalChainTopologicProperties = new ChainTopologicProperties(this,
                conventionalCalculations.averageGraphPathLength(),
                conventionalCalculations.clusteringCoefficient());

        chain.getFeatureContainer().addFeature(new ChainTopologicPropertiesContainer(this,
                fullPlipChainTopologicProperties,
                hydrogenPlipChainTopologicProperties,
                hydrophobicPlipChainTopologicProperties,
                conventionalChainTopologicProperties));

        chain.aminoAcids()
                .forEach(aminoAcid -> processInternally(aminoAcid,
                        fullPlipCalculations,
                        hydrogenPlipCalculations,
                        hydrophobicPlipCalculations,
                        conventionalCalculations));
    }

    private void processInternally(AminoAcid aminoAcid,
                                   ProteinGraphCalculations fullPlipCalculate,
                                   ProteinGraphCalculations hydrogenPlipCalculate,
                                   ProteinGraphCalculations hydrophobicPlipCalculate,
                                   ProteinGraphCalculations conventionalCalculate) {
        ResidueTopologicProperties fullPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                fullPlipCalculate.betweenness(aminoAcid),
                fullPlipCalculate.closeness(aminoAcid),
                fullPlipCalculate.clusteringCoefficient(aminoAcid),
                fullPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties hydrogenPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                hydrogenPlipCalculate.betweenness(aminoAcid),
                hydrogenPlipCalculate.closeness(aminoAcid),
                hydrogenPlipCalculate.clusteringCoefficient(aminoAcid),
                hydrogenPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties hydrophobicPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                hydrophobicPlipCalculate.betweenness(aminoAcid),
                hydrophobicPlipCalculate.closeness(aminoAcid),
                hydrophobicPlipCalculate.clusteringCoefficient(aminoAcid),
                hydrophobicPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties conventionalResidueTopologicProperties = new ResidueTopologicProperties(this,
                conventionalCalculate.betweenness(aminoAcid),
                conventionalCalculate.closeness(aminoAcid),
                conventionalCalculate.clusteringCoefficient(aminoAcid),
                conventionalCalculate.distinctNeighborhoodCount(aminoAcid));
        aminoAcid.getFeatureContainer().addFeature(new ResidueTopologicPropertiesContainer(this,
                fullPlipResidueTopologicProperties,
                hydrogenPlipResidueTopologicProperties,
                hydrophobicPlipResidueTopologicProperties,
                conventionalResidueTopologicProperties));
    }

    static class ProteinGraphCalculations {
        private final ProteinGraph graph;
        private final List<AminoAcid> nodes;
        private final Map<AminoAcid, ShortestPathAlgorithm.SingleSourcePaths<AminoAcid, DefaultEdge>> shortestPaths;
        private final DijkstraShortestPath<AminoAcid, DefaultEdge> dijkstraShortestPath;
        private final int numberOfNodes;
        private final double numberOfNodePairs;

        ProteinGraphCalculations(ProteinGraph graph) {
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
        double averageGraphPathLength() {
            double np = numberOfNodes * (numberOfNodes - 1) * 0.5;
            return 1 / np * graph.vertexSet()
                    .stream()
                    .mapToDouble(this::averageGraphPathLength)
                    .average()
                    .orElseThrow(() -> new IllegalArgumentException("could not compute average path length"));
        }

        double averageGraphPathLength(AminoAcid source) {
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
        double betweenness(AminoAcid node) {
            return  shortestPathsPassingThrough(node).size() / numberOfNodePairs;
        }

        /**
         * The maximal path length from this node to any other node within the graph.
         * See Amitai, 2004 for definition.
         * @param node the node to evaluate
         * @return the maximal path length involving this node
         */
        double closeness(AminoAcid node) {
            return 1 / graph.vertexSet().stream()
                    .map(target -> determineShortestPath(node, target))
                    .mapToInt(GraphPath::getLength)
                    .average()
                    .orElseThrow(() ->  new IllegalArgumentException("cannot evaluate closeness as graph is not fully connected"));
        }

        /**
         * Evaluates all nodes of a given graph by comparing the actual number of edges in the network of direct
         * neighbors to the theoretical maximum.
         * See Vendruscolo, 2002 for definition.
         * @return the clustering coefficient of this graph
         */
        double clusteringCoefficient() {
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
        double clusteringCoefficient(AminoAcid node) {
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

        int distinctNeighborhoodCount(AminoAcid aminoAcid) {
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
            return shortestPaths.get(source).getPath(target);
        }
    }
}
