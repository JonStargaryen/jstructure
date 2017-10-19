package de.bioforscher.jstructure.mathematics.graph.partitioning.algorithms;

import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;

import java.util.*;
import java.util.function.ToDoubleFunction;

/**
 * Implemented from Stijn van Dongen's (author of MCL algorithm) documentation: http://micans.org/mcl/
 * lecture notes: https://www.cs.ucsb.edu/~xyan/classes/CS595D-2009winter/MCL_Presentation2.pdf
 * and cytoscape's JS port: https://github.com/cytoscape/cytoscape.js-markov-cluster/blob/master/cytoscape-markov-cluster.js
 */
public class MCL implements GraphPartitioningAlgorithm {
    public static final double DEFAULT_EXPAND_FACTOR = 2;
    public static final double DEFAULT_INFLATE_FACTOR = 2;
    public static final double DEFAULT_MULT_FACTOR = 1;
    public static final double DEFAULT_MAX_ITERATIONS = 20;
    public static final ToDoubleFunction<Edge> DEFAULT_SIMILARITY_FUNCTION = Edge::getWeight;

    private final double expandFactor;
    private final double inflateFactor;
    private final double multFactor;
    private final double maxIterations;
    private final ToDoubleFunction<Edge> similarityFunction;

    public MCL(double expandFactor, double inflateFactor, double multFactor, double maxIterations, ToDoubleFunction<Edge> similarityFunction) {
        this.expandFactor = expandFactor;
        this.inflateFactor = inflateFactor;
        this.multFactor = multFactor;
        this.maxIterations = maxIterations;
        this.similarityFunction = similarityFunction;
    }

    public MCL() {
        this(DEFAULT_EXPAND_FACTOR,
                DEFAULT_INFLATE_FACTOR,
                DEFAULT_MULT_FACTOR,
                DEFAULT_MAX_ITERATIONS,
                DEFAULT_SIMILARITY_FUNCTION);
    }

    @Override
    public <N> PartitionedGraph<N> partitionGraph(Graph<N> graph) {
        List<N> nodes = graph.getNodes();
        List<Edge<N>> edges = graph.getEdges();

        // Map each node to its position in node array
        Map<N, Integer> id2Position = new HashMap<>();
        for (N node : nodes) {
            id2Position.put(node, id2Position.size());
        }

        // Generate stochastic matrix M from input graph G (should be symmetric/undirected)
        int n = nodes.size();
        int n2 = n * n;
        double[] M = new double[n2];
        double[] _M;

        for (Edge<N> edge : edges) {
            int i = id2Position.get(edge.getLeft());
            int j = id2Position.get(edge.getRight());

            double sim = getSimilarity(edge);

            // G should be symmetric and undirected
            M[i * n + j] += sim;
            M[j * n + i] += sim;
        }

        // Begin Markov cluster algorithm

        // Step 1: Add self loops to each node, ie. add multFactor to matrix diagonal
        addLoops(M, n);

        // Step 2: M = normalize( M );
        normalize(M, n);

        boolean isStillMoving = true;
        int iterations = 0;
        while(isStillMoving && iterations < maxIterations) {
            isStillMoving = false;

            // Step 3:
            _M = expand(M, n);

            // Step 4:
            M = inflate(_M, n);

            // Step 5: check to see if ~steady state has been reached
            if (!hasConverged(M, _M, n2, 4)) {
                isStillMoving = true;
            }

            iterations++;
        }

        return assign(M, n, graph);
    }

    private <N> PartitionedGraph<N> assign(double[] M, int n, Graph<N> graph) {
        List<Module<N>> clusters = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            List<N> cluster = new ArrayList<>();
            for (int j = 0; j < n; j++) {
                // Row-wise attractors and elements that they attract belong in same cluster
                if(Math.round(M[i * n + j] * 1000) / 1000 > 0) {
                    cluster.add(graph.getNodes().get(j));
                }
            }

            if(!cluster.isEmpty()) {
                clusters.add(new Module<>(String.valueOf(clusters.size() + 1), cluster));
            }
        }

        return new PartitionedGraph<>(graph, clusters);
    }

    private boolean hasConverged(double[] M, double[] _M, int n2, int roundFactor) {
        // Check that both matrices have the same elements (i,j)
        for (int i = 0; i < n2; i++ ) {
            // truncate to 'roundFactor' decimal places
            if (!numberEqual(M[i], _M[i], roundFactor)) {
                return false;
            }
        }
        return true;
    }

    boolean numberEqual(double v1, double v2, int precision) {
        v1 = Math.round(v1 * Math.pow(10, precision)) / Math.pow(10, precision);
        v2 = Math.round(v2 * Math.pow(10, precision)) / Math.pow(10, precision);
        return v1 == v2;
    }

    private double[] inflate(double[] M, int n) {
        double[] _M = new double[n * n];

        // M(i,j) ^ inflatePower
        for (int i = 0; i < n * n; i++) {
            _M[i] = Math.pow(M[i], inflateFactor);
        }

        normalize(_M, n);

        return _M;
    }

    private double[] expand(double[] M, int n) {
        double[] _M = Arrays.copyOf(M, M.length);

        for (int p = 0; p < expandFactor; p++) {
            M = mmult(M, _M, n);
        }

        return M;
    }

    private double[] mmult(double[] A, double[] B, int n) {
        double[] C = new double[n * n];

        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    C[i * n + j] += A[i * n + k] * B[k * n + j];
                }
            }
        }

        return C;
    }

    private void normalize(double[] M, int n) {
        double sum;
        for (int col = 0; col < n; col++) {
            sum = 0;
            for (int row = 0; row < n; row++) {
                sum += M[row * n + col];
            }
            for (int row = 0; row < n; row++) {
                M[row * n + col] = M[row * n + col] / sum;
            }
        }
    }

    private void addLoops(double[] M, int n) {
        for (int i = 0; i < n; i++) {
            M[i * n + i] = multFactor;
        }
    }

    private <N> double getSimilarity(Edge<N> edge) {
        return similarityFunction.applyAsDouble(edge);
    }
}
