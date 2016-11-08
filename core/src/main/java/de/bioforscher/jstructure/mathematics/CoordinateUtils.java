package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.function.Consumer;

/**
 * A collection of commonly used function to manipulate atom coordinates. In contrast to the the functions in
 * {@link LinearAlgebra3D} most of the function provided here are somewhat high-level and are not designed to deal with
 * <code>double[]</code> but rather with {@link Atom} objects. Also they will directly modify {@link Atom} instances
 * handed over to these functions.
 * @see LinearAlgebra3D
 * Created by S on 30.09.2016.
 */
public class CoordinateUtils {
    static final Logger logger = LoggerFactory.getLogger(CoordinateUtils.class);

    /**
     * Computes the center of mass of a collection of atoms. In contrast to the centroid, the center of mass considers
     * the unique mass of each atom in the structure (i.e. a {@link de.bioforscher.jstructure.model.structure.Element#H})
     * has a different mass than say {@link de.bioforscher.jstructure.model.structure.Element#C} and, thus, will affect
     * the center of mass computation in a less impactful way.
     * @param atomContainer a collection of atoms
     * @return the center of mass
     * @see CoordinateUtils#centroid(AtomContainer)
     */
    public static double[] centerOfMass(AtomContainer atomContainer) {
        return atomContainer.atoms()
                .collect(CenterOfMassAverager::new,
                         CenterOfMassAverager::accept,
                         CenterOfMassAverager::combine)
                .average();
    }

    static class CenterOfMassAverager implements Consumer<Atom> {
        private double[] coordinate = new double[3];
        private double mass = 0;

        double[] average() {
            return mass > 0 ? LinearAlgebra3D.divide(coordinate, mass) : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            coordinate = LinearAlgebra3D.add(coordinate, LinearAlgebra3D.multiply(atom.getCoordinates(), atom.getElement().getAtomicMass()));
            mass += atom.getElement().getAtomicMass();
        }

        void combine(CenterOfMassAverager other) {
            coordinate = LinearAlgebra3D.add(coordinate, other.coordinate);
            mass += other.mass;
        }
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atomContainer a collection of atoms
     * @return the maximal distance occurring between the centroid and any other atom
     * @see CoordinateUtils#maximalExtent(AtomContainer, double[])
     */
    public static double maximalExtent(AtomContainer atomContainer) {
        return maximalExtent(atomContainer, centroid(atomContainer));
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atomContainer a collection of atoms
     * @param centroid the centroid of this atom collection
     * @return the maximal distance occurring between the centroid and any other atom
     */
    public static double maximalExtent(AtomContainer atomContainer, final double[] centroid) {
        return atomContainer.atoms()
                .filter(AtomNameFilter.CA_ATOM_FILTER)
                .mapToDouble(atom -> LinearAlgebra3D.distance(atom.getCoordinates(), centroid))
                .max()
                .orElseThrow(() -> new IllegalArgumentException("cannot compute maximal extent for single atom"));
    }

    /**
     * Converts a collection of atoms to a <tt>n x 3</tt> matrix, where <tt>n</tt> is equal to the number of processed
     * atoms.
     * @param atomContainer a collection of atoms
     * @return a matrix containing the coordinates of all atoms
     */
    public static RealMatrix convertToMatrix(AtomContainer atomContainer) {
        List<Atom> atomList = atomContainer.getAtoms();
        double[][] matrix = new double[atomList.size()][3];
        for (int index = 0; index < atomList.size(); index++) {
            matrix[index] = atomList.get(index).getCoordinates();
        }
        return new Array2DRowRealMatrix(matrix);
    }

    /**
     * Employs a transformation (i.e. a rotation & translation or 'rototranslation') on a collection of atoms. The
     * operation will manipulate the provided atom's coordinates directly (rather than creating new Objects).
     * @param atomContainer a collection of atoms
     * @param translation the 3D vector describing the translation
     * @param rotation the <tt>3 x 3</tt> matrix describing the desired rotation
     */
    public static void transform(AtomContainer atomContainer, double[] translation, double[][] rotation) {
        atomContainer.atoms().forEach(new Transformation(translation, rotation)::transformCoordinates);
    }

    public static class Transformation {
        private double[] translation;
        private double[][] rotation;

        /**
         * Mere shift.
         * @param translation the translation vector
         */
        public Transformation(double[] translation) {
            this.translation = translation;
        }

        /**
         * Rototranslation.
         * @param translation the translation vector
         * @param rotation the rotation matrix
         */
        public Transformation(double[] translation, double[][] rotation) {
            this.translation = translation;
            this.rotation = rotation;
        }

        /**
         * Mere rotation.
         * @param rotation the rotation matrix
         */
        public Transformation(double[][] rotation) {
            this.rotation = rotation;
        }

        public Atom transformCoordinates(Atom atom) {
            double[] vector = atom.getCoordinates();
//            logger.debug("translated vector from {}", Arrays.toString(vector));
            // apply transformation if needed
            if(translation != null) {
                atom.setCoordinates(new double[] {
                    vector[0] + translation[0],
                    vector[1] + translation[1],
                    vector[2] + translation[2]
                });
            }

            vector = atom.getCoordinates();
//            logger.debug("\tto {}", Arrays.toString(vector));
            // apply rotation if needed
            if(rotation != null) {
                atom.setCoordinates(new double[] {
                    (rotation[0][0] * vector[0] + rotation[0][1] * vector[1] + rotation[0][2] * vector[2]),
                    (rotation[1][0] * vector[0] + rotation[1][1] * vector[1] + rotation[1][2] * vector[2]),
                    (rotation[2][0] * vector[0] + rotation[2][1] * vector[1] + rotation[2][2] * vector[2])
                });
            }
//            logger.debug("\tand rotated to {}", Arrays.toString(atom.getCoordinates()));
            return atom;
        }
    }

    /**
     * Centers a collection of atoms. This is achieved by computing the centroid of these points and subtracting this
     * point from each atom's coordinates. The operation will manipulate the provided atom's coordinates directly
     * (rather than creating new Objects).
     * @param atomContainer the collection of atoms to center
     * @see CoordinateUtils#centroid(AtomContainer)
     * @see CoordinateUtils#shift(AtomContainer, double[])
     */
    public static void center(AtomContainer atomContainer) {
        // invert the centroid/shift vector as it will be added to the coordinates by the shift function
        shift(atomContainer, LinearAlgebra3D.multiply(centroid(atomContainer), -1.0));
    }

    /**
     * Shifts (i.e. adding) each atom's coordinates by a defined shift vector. Multiply the shift vector by
     * <tt>-1.0</tt> to reverse direction. The operation will manipulate the provided atom's coordinates directly
     * (rather than creating new Objects).
     * @param atomContainer a collection of atoms to process
     * @param shiftVector the vector to add to each atom's coordinates
     */
    public static void shift(AtomContainer atomContainer, final double[] shiftVector) {
        atomContainer.atoms().forEach(new Transformation(shiftVector)::transformCoordinates);
    }

    /**
     * Computes the centroid (not center of mass!) of a collection of atoms.
     * @param atomContainer the atoms to process
     * @return a <code>double[]</code> describing the coordinates of the centroid
     */
    public static double[] centroid(AtomContainer atomContainer) {
        return atomContainer.atoms()
                .collect(CentroidAverager::new,
                         CentroidAverager::accept,
                         CentroidAverager::combine)
                .average();
    }

    static class CentroidAverager implements Consumer<Atom> {
        private double[] total = new double[3];
        private int count = 0;

        double[] average() {
            return count > 0 ? LinearAlgebra3D.divide(total, count) : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            total = LinearAlgebra3D.add(total, atom.getCoordinates());
            count++;
        }

        void combine(CentroidAverager other) {
            total = LinearAlgebra3D.add(total, other.total);
            count += other.count;
        }
    }

    /**
     * Computes the root-mean-square deviation between 2 atom sets. Both atom containers must have the same number of
     * atoms.
     * @param atomContainer1 a collection of atoms - it is supposed to share its ordering with the 2nd container
     * @param atomContainer2 another collection of atoms
     * @return the RMSD value of this alignment
     * @throws IllegalArgumentException is thrown when the number of atoms does no match
     */
    public static double calculateRMSD(AtomContainer atomContainer1, AtomContainer atomContainer2) throws IllegalArgumentException {
        double sd = Pair.sequentialPairsOf(atomContainer1.getAtoms(), atomContainer2.getAtoms())
            .mapToDouble(pair -> LinearAlgebra3D.distanceFast(pair.getFirst().getCoordinates(), pair.getSecond().getCoordinates()))
            .sum();
        return Math.sqrt(sd / atomContainer1.getAtoms().size());
    }
}