package de.bioforscher.jstructure.mathematics.mds;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

/**
 * The more or less (more less than more though) generic implementation of MDS.
 *
 * <p>This can be used to map data of some dimension (really high for gene expression levels or quite low for distance
 * maps) to some other dimension which can obviously different and will usually be 2 or 3. This can be exploited to
 * placed atom coordinates for proteins in the 3D space merely by knowing their inter-atomic distances.</p>
 * Created by S on 27.09.2016.
 */
public class MultiDimensionalScaling {
    /**
     * the default dimension of the space the MDS will embed into<br />
     * for sake of protein structure reconstruction, we aim at m=3
     */
    private static final int DEFAULT_TARGET_DIMENSION = 3;

    /**
     * the number of data points (determines the number of projected data points as well as the dimension of eigenvectors)
     */
    private int numberOfDataPoints;

    /**
     * the target dimension to which data points will be mapped (normally, 2 or
     * 3) - also known as m
     */
    private int targetDimension;

    /**
     * the distance map provided as input, describes an arbitrary (yet
     * symmetric) distance between an element i and j
     */
    private RealMatrix distanceMap;

    /**
     * the matrix A ('proximityMap' or D²) whereby each entry derived from the
     * input distances a_ij = -0.5*d_ij^2 thus, dim(distanceMap) =
     * dim(proximityMap)
     */
    private RealMatrix proximityMap;

    /**
     * the matrix B ('centeringMap') b_ij = a_ij - a_i* - a_j* + a_** a_i* - row
     * average a_j* - column average a_** - overall average
     */
    private RealMatrix centeringMap;

    /**
     * the m (target dimension) normalized eigenvectors which can subsequently be used to compose the embedded data points
     */
    private List<double[]> normalizedEigenvectors;

    /**
     * embedding in the space of dimension m
     */
    private List<double[]> embedding;

    /**
     * transforms this distance map into an embedded in the 3D space
     *
     * @param distanceMap
     *            a symmetric distance map of dim <tt>nxn</tt> - so <tt>n</tt>
     *            data points (e.g. atoms) are described by these distances
     * @return a list of size <tt>n</tt> and each entry has default dimensionality 3
     */
    public List<double[]> computeEmbedding(double[][] distanceMap) {
        // wrap the call, so the JAMA matrix object is not only way to use this class
        return computeEmbedding(distanceMap, DEFAULT_TARGET_DIMENSION);
    }

    public List<double[]> computeEmbedding(double[][] distanceMap, int targetDimension) {
        return computeEmbedding(new Array2DRowRealMatrix(distanceMap), targetDimension);
    }

    public List<double[]> computeEmbedding(RealMatrix distanceMap, int targetDimension) {
        numberOfDataPoints = distanceMap.getRowDimension();
        this.targetDimension = targetDimension;
        if (this.targetDimension > numberOfDataPoints) {
            throw new IllegalArgumentException("target dimension must not exceed number of data points");
        }

        this.distanceMap = distanceMap;

        proximityMap = computeSquaredProximityMap(this.distanceMap);

        centeringMap = computeConfiguration(proximityMap);

        normalizedEigenvectors = new ArrayList<>();
        embedding = new ArrayList<>();
        EigenDecomposition evd = new EigenDecomposition(centeringMap);
        // we are looking for the m biggest eigenvalues - they are at the last elements of the matrix
        RealMatrix eigenvectors = evd.getV();
        // Matrix eigenvalues = evd.getD();
        double[] eigenvalues = evd.getRealEigenvalues();
        Map<Integer, Double> eigenvalueMap = new HashMap<>();
        for (int eigenvalueIndex = 0; eigenvalueIndex < eigenvalues.length; eigenvalueIndex++) {
            eigenvalueMap.put(eigenvalueIndex, eigenvalues[eigenvalueIndex]);
        }
        List<Entry<Integer, Double>> sortedEigenvalues = entriesSortedByValues(eigenvalueMap).subList(0, targetDimension);

        // normalize eigenvectors
        for (Entry<Integer, Double> sortedEigenvalue : sortedEigenvalues) {
            if (sortedEigenvalue.getValue() <= 0) {
                throw new IllegalArgumentException("eigenvalue is negative: " + sortedEigenvalue.getValue());
            }
            double[] eigenvector = eigenvectors.getColumn(sortedEigenvalue.getKey());
            normalizedEigenvectors.add(normalize(eigenvector, Math.sqrt(sortedEigenvalue.getValue())));
        }

        // compose embedded data points from normalized eigenvectors
        for (int dataPointIndex = 0; dataPointIndex < numberOfDataPoints; dataPointIndex++) {
            double[] dataPoint = new double[this.targetDimension];
            for (int dataPointDimension = 0; dataPointDimension < this.targetDimension; dataPointDimension++) {
                dataPoint[dataPointDimension] = normalizedEigenvectors.get(dataPointDimension)[dataPointIndex];
            }
            embedding.add(dataPoint);
        }

        return embedding;
    }

    // TODO: maybe move to some Utils class
    private static <K, V extends Comparable<? super V>> List<Entry<K, V>> entriesSortedByValues(Map<K, V> map) {
        List<Entry<K, V>> sortedEntries = new ArrayList<>(map.entrySet());
        sortedEntries.sort((o1, o2) -> o2.getValue().compareTo(o1.getValue()));
        return sortedEntries;
    }

    private double[] normalize(double[] vector, double normalizationFactor) {
        for (int dimension = 0; dimension < vector.length; dimension++) {
            vector[dimension] *= normalizationFactor;
        }
        return vector;
    }

    private RealMatrix computeConfiguration(RealMatrix proximityMap) {
        RealMatrix centeringMap = new Array2DRowRealMatrix(numberOfDataPoints, numberOfDataPoints);

        double[] rowAverage = new double[numberOfDataPoints];
        double[] columnAverage = new double[numberOfDataPoints];
        double overallAverage = 0;
        // assess rows and overall average
        for (int row = 0; row < numberOfDataPoints; row++) {
            double tempRowAverage = 0;
            for (int column = 0; column < numberOfDataPoints; column++) {
                double entry = proximityMap.getEntry(row, column);
                tempRowAverage += entry;
                overallAverage += entry;
            }
            rowAverage[row] = tempRowAverage / numberOfDataPoints;
        }
        overallAverage /= numberOfDataPoints * numberOfDataPoints;

        // assess columns
        for (int column = 0; column < numberOfDataPoints; column++) {
            double tempColumnAverage = 0;
            for (int row = 0; row < numberOfDataPoints; row++) {
                tempColumnAverage += proximityMap.getEntry(row, column);
            }
            columnAverage[column] = tempColumnAverage / numberOfDataPoints;
        }

        for (int row = 0; row < numberOfDataPoints; row++) {
            for (int column = 0; column < numberOfDataPoints; column++) {
                // b_ij = a_ij - a_i* - a_j* + a_**
                centeringMap.setEntry(row, column,
                        proximityMap.getEntry(row, column) - rowAverage[row] - columnAverage[column] + overallAverage);
            }
        }
        return centeringMap;
    }

    private RealMatrix computeSquaredProximityMap(RealMatrix distanceMap) {
        RealMatrix proximityMap = new Array2DRowRealMatrix(numberOfDataPoints, numberOfDataPoints);

        for (int row = 0; row < numberOfDataPoints; row++) {
            for (int column = 0; column < numberOfDataPoints; column++) {
                // a_ij = -0.5*d_ij^2
                double entry = distanceMap.getEntry(row, column);
                proximityMap.setEntry(row, column, -0.5 * entry * entry);
            }
        }
        return proximityMap;
    }
}