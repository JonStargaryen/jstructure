package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.concurrent.ThreadSafe;
import java.util.List;

import static de.bioforscher.jstructure.StandardFormat.format;

/**
 * This class provides access to basic linear algebraic operations on either <code>double[]</code> or model instances
 * providing fields with e.g. coordinates. Instances are immutable.
 * Created by bittrich on 4/28/17.
 */
@ThreadSafe
public class LinearAlgebra {
    private static final Logger logger = LoggerFactory.getLogger(LinearAlgebra.class);

    private LinearAlgebra() {
        // deny direct instantiation
    }

    /**
     * Construct a wrapper for basic operations on a 3D double array.
     * @param value the value to perform algebraic operations on
     * @return an wrapping instance for algebraic operations
     */
    public static PrimitiveDoubleArrayLinearAlgebra on(double[] value) {
        return new PrimitiveDoubleArrayLinearAlgebra(value);
    }

    /**
     * Construct a wrapper for basic operations on an atom instance.
     * @param atom the atom to perform algebraic operations on
     * @return an wrapping instance for algebraic operations
     */
    public static PrimitiveDoubleArrayLinearAlgebra on(Atom atom) {
        return new PrimitiveDoubleArrayLinearAlgebra(atom.getCoordinates());
    }

    /**
     * Construct a wrapper for basic operations on an atom container.
     * @param atomContainer the atom to perform algebraic operations on
     * @return an wrapping instance for algebraic operations
     */
    public static AtomContainerLinearAlgebra on(AtomContainer atomContainer) {
        return new AtomContainerLinearAlgebra(atomContainer);
    }

    public static <N> GraphLinearAlgebra<N> on(Graph<N> graph) {
        return new GraphLinearAlgebra<>(graph);
    }

    /**
     * Caps a number to a defined interval.
     * @param lowerBound the lowest accepted value x_l
     * @param value the raw value
     * @param upperBound the highest accepted value x_u
     * @return the original value if it lies in the interval [x_l,x_u], else x_l or x_u
     */
    static double capToInterval(double lowerBound, double value, double upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }

    /**
     * Caps a number to a defined interval.
     * @param lowerBound the lowest accepted value x_l
     * @param value the raw value
     * @param upperBound the highest accepted value x_u
     * @return the original value if it lies in the interval [x_l,x_u], else x_l or x_u
     */
    static int capToInterval(int lowerBound, int value, int upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }

    public interface AlgebraicOperations {

    }

    public static class PrimitiveDoubleArrayLinearAlgebra implements AlgebraicOperations {
        private final double x;
        private final double y;
        private final double z;

        PrimitiveDoubleArrayLinearAlgebra(double... value) {
            this.x = value[0];
            this.y = value[1];
            this.z = value[2];
        }

        public PrimitiveDoubleArrayLinearAlgebra add(Atom vectorToAdd) {
            return add(vectorToAdd.getCoordinates());
        }

        public PrimitiveDoubleArrayLinearAlgebra add(double[] vectorToAdd) {
            return new PrimitiveDoubleArrayLinearAlgebra(x + vectorToAdd[0],
                    y + vectorToAdd[1],
                    z + vectorToAdd[2]);
        }

        public PrimitiveDoubleArrayLinearAlgebra add(PrimitiveDoubleArrayLinearAlgebra vectorToAdd) {
            return add(vectorToAdd.getValue());
        }

        public double angle(Atom secondVector) {
            return angle(secondVector.getCoordinates());
        }

        public double angle(double[] secondVector) {
            double vDot = dotProduct(secondVector) / (norm() * LinearAlgebra.on(secondVector).norm());
            return Math.toDegrees(Math.acos(capToInterval(-1, vDot, 1)));
        }

        public double angle(PrimitiveDoubleArrayLinearAlgebra secondVector) {
            return angle(secondVector.getValue());
        }

        public double distance(Atom secondVector) {
            return distance(secondVector.getCoordinates());
        }

        public double distance(double[] secondVector) {
            return Math.sqrt(distanceFast(secondVector));
        }

        public double distance(PrimitiveDoubleArrayLinearAlgebra secondVector) {
            return distance(secondVector.getValue());
        }

        public double distanceFast(Atom secondVector) {
            return distanceFast(secondVector.getCoordinates());
        }

        public double distanceFast(double[] secondVector) {
            return (x - secondVector[0]) * (x - secondVector[0]) +
                    (y - secondVector[1]) * (y - secondVector[1]) +
                    (z - secondVector[2]) * (z - secondVector[2]);
        }

        public double distanceFast(PrimitiveDoubleArrayLinearAlgebra secondVector) {
            return distanceFast(secondVector.getValue());
        }

        public PrimitiveDoubleArrayLinearAlgebra divide(double divisor) {
            return new PrimitiveDoubleArrayLinearAlgebra(x / divisor,
                    y / divisor,
                    z / divisor);
        }

        public double dotProduct(Atom secondVector) {
            return dotProduct(secondVector.getCoordinates());
        }

        public double dotProduct(double[] secondVector) {
            return x * secondVector[0] +
                    y * secondVector[1] +
                    z * secondVector[2];
        }

        public double dotProduct(PrimitiveDoubleArrayLinearAlgebra secondVector) {
            return dotProduct(secondVector.getValue());
        }

        public double[] getValue() {
            return new double[] { x, y, z };
        }

        public PrimitiveDoubleArrayLinearAlgebra multiply(double scalar) {
            return new PrimitiveDoubleArrayLinearAlgebra(x * scalar,
                    y * scalar,
                    z * scalar);
        }

        public PrimitiveDoubleArrayLinearAlgebra multiply(double[][] matrix) {
            double tx = matrix[0][0] * x + matrix[1][0] * y + matrix[2][0] * z;
            double ty = matrix[0][1] * x + matrix[1][1] * y + matrix[2][1] * z;
            double tz = matrix[0][2] * x + matrix[1][2] * y + matrix[2][2] * z;
            return new PrimitiveDoubleArrayLinearAlgebra(tx,
                    ty,
                    tz);
        }

        public double norm() {
            return Math.sqrt(dotProduct(getValue()));
        }

        public PrimitiveDoubleArrayLinearAlgebra normalize() {
            return divide(norm());
        }

        public PrimitiveDoubleArrayLinearAlgebra subtract(Atom vectorToSubtract) {
            return subtract(vectorToSubtract.getCoordinates());
        }

        public PrimitiveDoubleArrayLinearAlgebra subtract(double[] vectorToSubtract) {
            return new PrimitiveDoubleArrayLinearAlgebra(x - vectorToSubtract[0],
                    y - vectorToSubtract[1],
                    z - vectorToSubtract[2]);
        }

        public PrimitiveDoubleArrayLinearAlgebra subtract(PrimitiveDoubleArrayLinearAlgebra vectorToSubtract) {
            return subtract(vectorToSubtract.getValue());
        }

        public PrimitiveDoubleArrayLinearAlgebra vectorProduct(Atom secondVector) {
            return vectorProduct(secondVector.getCoordinates());
        }

        public PrimitiveDoubleArrayLinearAlgebra vectorProduct(double[] secondVector) {
            double tx = y * secondVector[2] - z * secondVector[1];
            double ty = z * secondVector[0] - x * secondVector[2];
            double tz = x * secondVector[1] - y * secondVector[0];
            return new PrimitiveDoubleArrayLinearAlgebra(tx,
                    ty,
                    tz);
        }

        public PrimitiveDoubleArrayLinearAlgebra vectorProduct(PrimitiveDoubleArrayLinearAlgebra secondVector) {
            return vectorProduct(secondVector.getValue());
        }

        @Override
        public String toString() {
            return "[" + format(x) + ", " + format(y) + ", " + format(z) + "]";
        }
    }

    public static class AtomContainerLinearAlgebra implements AlgebraicOperations {
        private final AtomContainer atomContainer;

        AtomContainerLinearAlgebra(AtomContainer atomContainer) {
            this.atomContainer = atomContainer;
        }

        /**
         * Centers a collection of atomContainer. This is achieved by computing the centroid of these points and subtracting this
         * point from each atom's coordinates. The operation will manipulate the provided atom's coordinates directly
         * (rather than creating new Objects).
         * @see AtomContainerLinearAlgebra#centroid()
         * @return provides the centroid of this entity
         */
        public PrimitiveDoubleArrayLinearAlgebra center() {
            PrimitiveDoubleArrayLinearAlgebra centroid = centroid();
            // invert the centroid/shift vector as it will be added to the coordinates by the shift function
            transform(centroid().multiply(-1).getValue());
            return centroid;
        }

        /**
         * Computes the center of mass of a collection of atomContainer. In contrast to the centroid, the center of mass considers
         * the unique mass of each atom in the structure (i.e. a {@link de.bioforscher.jstructure.model.structure.Element#H})
         * has a different mass than say {@link de.bioforscher.jstructure.model.structure.Element#C} and, thus, will affect
         * the center of mass computation in a less impacting way.
         * @return the center of mass
         * @see AtomContainerLinearAlgebra#centroid()
         */
        public PrimitiveDoubleArrayLinearAlgebra centerOfMass() {
            return new PrimitiveDoubleArrayLinearAlgebra(atomContainer.atoms().collect(StructureCollectors.toCenterOfMass()));
        }

        /**
         * Computes the centroid (not center of mass!) of a collection of atomContainer.
         * @return the coordinates of the centroid
         */
        public PrimitiveDoubleArrayLinearAlgebra centroid() {
            return new PrimitiveDoubleArrayLinearAlgebra(atomContainer.atoms().collect(StructureCollectors.toCentroid()));
        }

        /**
         * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
         * @return the maximal distance occurring between the centroid and any other atom
         */
        public double maximalExtent() {
            return maximalExtent(centroid().getValue());
        }

        /**
         * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
         * @param centroid the centroid of this atom collection
         * @return the maximal distance occurring between the centroid and any other atom
         */
        public double maximalExtent(final double[] centroid) {
            return Math.sqrt(atomContainer.atoms()
                    .mapToDouble(atom -> LinearAlgebra.on(atom.getCoordinates()).distanceFast(centroid))
                    .max()
                    .orElseThrow(() -> new IllegalArgumentException("cannot calculate maximal extent for single atom")));
        }

        public AtomContainerLinearAlgebra transform(double[] translation) {
            return transform(new Transformation(translation));
        }

        public AtomContainerLinearAlgebra transform(double[][] rotation) {
            return transform(new Transformation(rotation));
        }

        public AtomContainerLinearAlgebra transform(double[] translation, double[][] rotation) {
            return transform(new Transformation(translation, rotation));
        }

        public AtomContainerLinearAlgebra transform(Transformation transformation) {
            atomContainer.atoms().forEach(transformation::transform);
            return this;
        }
    }

    public static class GraphLinearAlgebra<N> implements AlgebraicOperations {
        private final Graph<N> graph;
        private final double numberOfNodePairs;
        private int[][] dist;

        GraphLinearAlgebra(Graph<N> graph) {
            this.graph = graph;
            this.numberOfNodePairs = graph.getNumberOfNodes() * (graph.getNumberOfNodes() - 1) * 0.5;
            this.dist = new int[graph.getNumberOfNodes()][graph.getNumberOfNodes()];

            int i, j, k;

            for (i = 0; i < graph.getNumberOfNodes(); i++) {
                for (j = 0; j < graph.getNumberOfNodes(); j++) {
                    if(i == j) {
                        dist[i][j] = 0;
                    } else {
                        //TODO support for weighted graphs
                        dist[i][j] = graph.containsEdge(graph.getNodes().get(i), graph.getNodes().get(j)) ?
                                1 : 999999;
                    }
                }
            }

            for (k = 0; k < graph.getNumberOfNodes(); k++) {
                for (i = 0; i < graph.getNumberOfNodes(); i++)  {
                    for (j = 0; j < graph.getNumberOfNodes(); j++) {
                        // If vertex k is on the shortest path from
                        // i to j, then update the value of dist[i][j]
                        if (dist[i][k] + dist[k][j] < dist[i][j])
                            dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }

        /**
         * Determines the average path length between all pairs of nodes in a given graph. Graph must be connected.
         * See Vendruscolo, 2002 for definition.
         * @return the average path length of this graph
         */
        public double averageGraphPathLength() {
            int numberOfNodes = graph.getNumberOfNodes();
            double np = numberOfNodes * (numberOfNodes - 1) * 0.5;
            return 1 / np * SetOperations.uniquePairsOf(graph.getNodes())
                    .mapToInt(this::determineShortestPathLength)
                    .sum();
        }

        /**
         * Betweenness is defined as the number shortest paths on the graph passing through this node normalized by the
         * total number of pairs of nodes of the graph.
         * See Vendruscolo, 2002 for definition.
         * @param node the node to evaluate
         * @return the betweenness of this node
         */
        public double betweenness(N node) {
            return 0;
//            return shortestPathsPassingThrough(node) / numberOfNodePairs;
        }

        /**
         * The maximal path length from this node to any other node within the graph.
         * See Amitai, 2004 for definition.
         * @param node the node to evaluate
         * @return the maximal path length involving this node
         */
        public int closeness(N node) {
            return graph.nodes()
                    .filter(n -> !node.equals(n))
                    .map(n -> new Pair<>(node, n))
                    .mapToInt(this::determineShortestPathLength)
                    .max()
                    .orElseThrow(() ->  new IllegalArgumentException("cannot evaluate closeness as graph is not fully connected"));
        }

        private int determineShortestPathLength(Pair<N, N> pair) {
            int i = graph.getNodes().indexOf(pair.getLeft());
            int j = graph.getNodes().indexOf(pair.getRight());
            return dist[i][j];
        }

//        /**
//         * Determine the number of shortest paths present in the graph passing through the specified node.
//         * @param node the node to evaluate
//         * @return the number of shortest paths on the graph passing through
//         */
//        private int shortestPathsPassingThrough(N node) {
//            return (int) SetOperations.uniquePairsOf(graph.getNodes())
//                    .map(this::determineShortestPath)
//                    .filter(graphPath -> graphPath.contains(node))
//                    .count();
//        }

//        private LinkedList<N> determineShortestPath(Pair<N, N> pair) {
//            DijkstraAlgorithm dijkstraAlgorithm = new DijkstraAlgorithm(graph);
//            dijkstraAlgorithm.execute(pair.getLeft());
//            return dijkstraAlgorithm.getPath(pair.getRight());
//        }
//
//        public class DijkstraAlgorithm {
//            private final List<N> nodes;
//            private final List<Edge<N>> edges;
//            private Set<N> settledNodes;
//            private Set<N> unSettledNodes;
//            private Map<N, N> predecessors;
//            private Map<N, Double> distance;
//
//            public DijkstraAlgorithm(Graph<N> graph) {
//                // create a copy of the array so that we can operate on this array
//                this.nodes = new ArrayList<>(graph.getNodes());
//                this.edges = new ArrayList<>(graph.getEdges());
//            }
//
//            public void execute(N source) {
//                settledNodes = new HashSet<>();
//                unSettledNodes = new HashSet<>();
//                distance = new HashMap<>();
//                predecessors = new HashMap<>();
//                distance.put(source, 0.0);
//                unSettledNodes.add(source);
//                while (unSettledNodes.size() > 0) {
//                    N node = getMinimum(unSettledNodes);
//                    settledNodes.add(node);
//                    unSettledNodes.remove(node);
//                    findMinimalDistances(node);
//                }
//            }
//
//            private void findMinimalDistances(N node) {
//                List<N> adjacentNodes = graph.getNeighborsFor(node);
//                for(N target : adjacentNodes) {
//                    if (getShortestDistance(target) > getShortestDistance(node)
//                            + getDistance(node, target)) {
//                        distance.put(target, getShortestDistance(node) + getDistance(node, target));
//                        predecessors.put(target, node);
//                        unSettledNodes.add(target);
//                    }
//                }
//
//            }
//
//            private double getDistance(N node, N target) {
//                for (Edge edge : edges) {
//                    if (edge.contains(node) && edge.contains(target)) {
//                        return edge.getWeight();
//                    }
//                }
//                throw new RuntimeException("Should not happen");
//            }
//
//            private N getMinimum(Set<N> vertexes) {
//                N minimum = null;
//                for (N vertex : vertexes) {
//                    if (minimum == null) {
//                        minimum = vertex;
//                    } else {
//                        if (getShortestDistance(vertex) < getShortestDistance(minimum)) {
//                            minimum = vertex;
//                        }
//                    }
//                }
//                return minimum;
//            }
//
//            private boolean isSettled(N vertex) {
//                return settledNodes.contains(vertex);
//            }
//
//            private double getShortestDistance(N destination) {
//                Double d = distance.get(destination);
//                if (d == null) {
//                    return Integer.MAX_VALUE;
//                } else {
//                    return d;
//                }
//            }
//
//            /*
//             * This method returns the path from the source to the selected target and
//             * NULL if no path exists
//             */
//            public LinkedList<N> getPath(N target) {
//                LinkedList<N> path = new LinkedList<>();
//                N step = target;
//                // check if a path exists
//                if (predecessors.get(step) == null) {
//                    return null;
//                }
//                path.add(step);
//                while (predecessors.get(step) != null) {
//                    step = predecessors.get(step);
//                    path.add(step);
//                }
//                // Put it into the correct order
//                Collections.reverse(path);
//                return path;
//            }
//        }

//        private GraphPath<N> determineShortestPath(Pair<N, N> pair) {
//            return findBasedOnPredicate(pair.getLeft(), n -> n.equals(pair.getRight()));
//        }
//
//        /**
//         * Returns the shortest path originating from the source node. The target node is the first node that satisfies
//         * target predicate. E.g. To to search for a specific node in the graph it is possible to use the identifier in the
//         * predicate. If no path can be found null is returned.
//         *
//         * @param sourceNode The source node.
//         * @param targetPredicate The predicate the target has to fulfill.
//         * @return The shortest path.
//         */
//        public GraphPath<N> findBasedOnPredicate(N sourceNode, Predicate<N> targetPredicate) {
//            this.queue.clear();
//            this.distances.clear();
//            this.predecessors.clear();
//            this.queue.offer(sourceNode);
//            this.distances.put(sourceNode, 0);
//            // processes
//            while (!queue.isEmpty()) {
//                N currentNode = queue.poll();
//                for (N neighbour : graph.getNeighborsFor(currentNode)) {
//                    GraphPath<N> path = checkTarget(currentNode, neighbour, targetPredicate);
//                    if (path != null) {
//                        return path;
//                    }
//                }
//            }
//            return null;
//        }
//
//        /**
//         * Checks whether the current target fulfills the predicate. If this is the case the path from the source to the
//         * target is returned. Otherwise the target is added to the queue, it is referenced in the predecessors map and the
//         * distance is set in the distance map.
//         *
//         * @param source The source node.
//         * @param target The target node.
//         * @param targetPredicate The predicate to fulfill.
//         * @return The shortest path if the target fulfills the predicate and null otherwise.
//         */
//        private GraphPath<N> checkTarget(N source, N target, Predicate<N> targetPredicate) {
//            System.out.println("checking " + source + " and " + target);
//            if (!distances.containsKey(target)) {
//                // until predicate is fulfilled the first time
//                if (targetPredicate.test(target)) {
//                    predecessors.put(target, source);
//                    return getPath(target);
//                }
//                // calculate distance and offer to queue
//                distances.put(target, distances.get(source) + 1);
//                predecessors.put(target, source);
//                queue.offer(target);
//                System.out.println("processed " + source + ", put " + target + " on queue");
//            }
//            return null;
//        }
//
//        /**
//         * Builds the path to the given target node.
//         *
//         * @param targetNode The target node.
//         * @return The path to the node.
//         */
//        private GraphPath<N> getPath(N targetNode) {
//            List<N> path = new ArrayList<>();
//            N step = targetNode;
//            if (predecessors.get(step) == null) {
//                return null;
//            }
//            path.add(step);
//            while (predecessors.get(step) != null) {
//                step = predecessors.get(step);
//                path.add(step);
//            }
//            Collections.reverse(path);
//            return new GraphPath<>(path);
//        }
//
//        private int determineShortestPathLength(Pair<N, N> pair) {
//            LinkedList<N> graphPath = determineShortestPath(pair);
//            if(graphPath == null) {
//                return Integer.MAX_VALUE;
//            } else {
//                return graphPath.size();
//            }
//        }

        /**
         * Evaluates all nodes of a given graph by comparing the actual number of edges in the network of direct
         * neighbors to the theoretical maximum.
         * See Vendruscolo, 2002 for definition.
         * @return the clustering coefficient of this graph
         */
        public double clusteringCoefficient() {
            return graph.nodes()
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
        public double clusteringCoefficient(N node) {
            List<N> neighbors = graph.getNeighborsFor(node);
            int numberOfNeighbors = neighbors.size();
            if(numberOfNeighbors == 1) {
                return 0;
            }
            int actualNumberOfEdges = (int) SetOperations.uniquePairsOf(neighbors)
                    .filter(graph::containsEdge)
                    .count();
            return actualNumberOfEdges / (0.5 * (numberOfNeighbors * (numberOfNeighbors - 1)));
        }
    }
}
