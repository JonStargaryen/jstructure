package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.GraphPath;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import javax.annotation.concurrent.ThreadSafe;
import java.util.*;

import static de.bioforscher.jstructure.StandardFormat.format;

/**
 * This class provides access to basic linear algebraic operations on either <code>double[]</code> or model instances
 * providing fields with e.g. coordinates. Instances are immutable.
 * Created by bittrich on 4/28/17.
 */
@ThreadSafe
public class LinearAlgebra {
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

        GraphLinearAlgebra(Graph<N> graph) {
            this.graph = graph;
            this.numberOfNodePairs = graph.getNumberOfNodes() * (graph.getNumberOfNodes() - 1) * 0.5;
        }

        /**
         * Determines the average path length between all pairs of nodes in a given graph. Graph must be connected.
         * See Vendruscolo, 2002 for definition.
         * @return the average path length of this graph
         */
        public double graphPathLength() {
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
            return shortestPathsPassingThrough(node) / numberOfNodePairs;
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
                    .map(n -> new Pair<>(n, node))
                    .mapToInt(this::determineShortestPathLength)
                    .max()
                    .orElseThrow(() ->  new IllegalArgumentException("cannot evaluate closeness as graph is not fully connected"));
        }

        /**
         * Determine the number of shortest paths present in the graph passing through the specified node.
         * @param node the node to evaluate
         * @return the number of shortest paths on the graph passing through
         */
        private int shortestPathsPassingThrough(N node) {
            return (int) SetOperations.uniquePairsOf(graph.getNodes())
                    .map(this::determineShortestPath)
                    .filter(graphPath -> graphPath.getElements().contains(node))
                    .count();
        }

        private GraphPath<N> determineShortestPath(Pair<N, N> pair) {
            N source = pair.getLeft();
            N target = pair.getRight();

            // initialize maps and queue
            Map<N, Integer> distanceMap = new HashMap<>();
            Map<N, N> predecessors = new HashMap<>();
            Queue<N> queue = new LinkedList<>();
            queue.offer(source);
            distanceMap.put(source, 0);

            while(!queue.isEmpty()) {
                N node = queue.poll();
                for(N neighbor : graph.getNeighborsFor(node)) {
                    if(!distanceMap.containsKey(target)) {
                        predecessors.put(target, source);

                        if(neighbor.equals(target)) {
                            return reconstructShortestPath(predecessors, target);
                        }

                        distanceMap.put(neighbor, distanceMap.get(source) + 1);
                        queue.offer(neighbor);
                    }
                }
            }

            // no valid path found
            return new GraphPath<>(new ArrayList<>());
        }

        private GraphPath<N> reconstructShortestPath(Map<N, N> predecessors, N target) {
            N step = target;
            List<N> steps = new ArrayList<>();

            // return empty path
            if(!predecessors.containsKey(step)) {
                return new GraphPath<>(steps);
            }

            steps.add(step);
            while(predecessors.containsKey(step)) {
                step = predecessors.get(step);
                steps.add(step);
            }

            Collections.reverse(steps);
            return new GraphPath<>(steps);
        }

        private int determineShortestPathLength(Pair<N, N> pair) {
            return determineShortestPath(pair).size();
        }

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
