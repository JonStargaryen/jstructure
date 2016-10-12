package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.model.structure.filter.AtomPairDistanceCutoffFilter;
import de.bioforscher.jstructure.model.structure.filter.GroupPairDistanceCutoffFilter;
import de.bioforscher.jstructure.model.structure.filter.ResiduePairDistanceCutoffFilter;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A collection of commonly used function to manipulate atom coordinates. In contrast to the the functions in
 * {@link LinearAlgebra3D} most of the function provided here are somewhat high-level and are not designed to deal with
 * <code>double[]</code> but rather with {@link Atom} objects. Also they will directly modify {@link Atom} instances
 * handed over to these functions.
 * @see LinearAlgebra3D
 * Created by S on 30.09.2016.
 */
public class CoordinateUtils {
    /**
     * Returns a stream of all atom combination whose euclidean distance is below a defined threshold.
     * @param atoms a collection of atoms
     * @param distanceCutoff the distance threshold below two atoms are considered to be in contact
     * @return all atom pairs in contact according to the distance cutoff
     * @see Pair#unorderedPairsOf(List)
     * @see AtomPairDistanceCutoffFilter
     */
    public static Stream<Pair<Atom>> atomPairsInContact(Stream<Atom> atoms, double distanceCutoff) {
        final AtomPairDistanceCutoffFilter distanceCutoffFilter = new AtomPairDistanceCutoffFilter(distanceCutoff);
        return Pair.unorderedPairsOf(atoms.collect(Collectors.toList())).filter(distanceCutoffFilter);
    }

    /**
     * Computes the center of mass of a collection of atoms. In contrast to the centroid, the center of mass considers
     * the unique mass of each atom in the structure (i.e. a {@link de.bioforscher.jstructure.model.structure.Element#H})
     * has a different mass than say {@link de.bioforscher.jstructure.model.structure.Element#C} and, thus, will affect
     * the center of mass computation in a less impactful way.
     * @param atoms a collection of atoms
     * @return the center of mass
     * @see CoordinateUtils#centroid(Stream)
     */
    public static double[] centerOfMass(Stream<Atom> atoms) {
        //TODO implement
        throw new UnsupportedOperationException("not implemented");
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atoms a collection of atoms
     * @return the maximal distance occurring between the centroid and any other atom
     * @see CoordinateUtils#maximalExtent(Stream, double[])
     */
    public static double maximalExtent(Stream<Atom> atoms) {
        List<Atom> atomList = atoms.collect(Collectors.toList());
        return maximalExtent(atomList.stream(), centroid(atomList.stream()));
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atoms a collection of atoms
     * @param centroid the centroid of this atom collection
     * @return the maximal distance occurring between the centroid and any other atom
     */
    public static double maximalExtent(Stream<Atom> atoms, double[] centroid) {
        return atoms.filter(AtomNameFilter.CA_ATOM_FILTER).mapToDouble(atom ->
                LinearAlgebra3D.distance(atom.getCoordinates(), centroid)).max().getAsDouble();
    }

    /**
     * Converts a collection of atoms to a <tt>n x 3</tt> matrix, where <tt>n</tt> is equal to the number of processed
     * atoms.
     * @param atoms a collection of atoms
     * @return a matrix containing the coordinates of all atoms
     */
    public static RealMatrix convertToMatrix(Stream<Atom> atoms) {
        List<Atom> atomList = atoms.collect(Collectors.toList());
        double[][] matrix = new double[atomList.size()][3];
        for (int index = 0; index < atomList.size(); index++) {
            matrix[index] = atomList.get(index).getCoordinates();
        }
        return new Array2DRowRealMatrix(matrix);
    }

    /**
     * Employs a transformation (i.e. a rotation & translation or 'rototranslation') on a collection of atoms. The
     * operation will manipulate the provided atom's coordinates directly (rather than creating new Objects).
     * @param atoms a collection of atoms
     * @param translation the 3D vector describing the translation
     * @param rotation the <tt>3 x 3</tt> matrix describing the desired rotation
     */
    public static Stream<Atom> transform(Stream<Atom> atoms, double[] translation, double[][] rotation) {
        return atoms.peek(atom -> {
            double[] vector = atom.getCoordinates();
            atom.setCoordinates(new double[] {
                    (rotation[0][0] * vector[0] + rotation[0][1] * vector[1] + rotation[0][2] * vector[2]) + translation[0],
                    (rotation[1][0] * vector[0] + rotation[1][1] * vector[1] + rotation[1][2] * vector[2]) + translation[1],
                    (rotation[2][0] * vector[0] + rotation[2][1] * vector[1] + rotation[2][2] * vector[2]) + translation[2]
            });
        });
    }

    /**
     * Centers a collection of atoms. This is achieved by computing the centroid of these points and subtracting this
     * point from each atom's coordinates. The operation will manipulate the provided atom's coordinates directly
     * (rather than creating new Objects).
     * @param atoms the collection of atoms to center
     * @see CoordinateUtils#centroid(Stream)
     * @see CoordinateUtils#shift(Stream, double[])
     */
    public static Stream<Atom> center(Stream<Atom> atoms) {
        List<Atom> atomList = atoms.collect(Collectors.toList());
        // invert the centroid/shift vector as it will be added to the coordinates by the shift function
        return shift(atomList.stream(), LinearAlgebra3D.multiply(centroid(atomList.stream()), -1.0));
    }

    /**
     * Shifts (i.e. adding) each atom's coordinates by a defined shift vector. Multiply the shift vector by
     * <tt>-1.0</tt> to reverse direction. The operation will manipulate the provided atom's coordinates directly
     * (rather than creating new Objects).
     * @param atoms a collection of atoms to process
     * @param shiftVector the vector to add to each atom's coordinates
     */
    public static Stream<Atom> shift(Stream<Atom> atoms, double[] shiftVector) {
        return atoms.peek(atom -> {
            atom.setCoordinates(LinearAlgebra3D.add(atom.getCoordinates(), shiftVector));
        });
    }

    /**
     * Computes the centroid (not center of mass!) of a collection of atoms.
     * @param atoms the atoms to process
     * @return a <code>double[]</code> describing the coordinates of the centroid
     */
    public static double[] centroid(Stream<Atom> atoms) {
        List<Atom> atomList = atoms.collect(Collectors.toList());
        return new double[] {
                atomList.stream().mapToDouble(atom -> atom.getCoordinates()[0]).average().getAsDouble(),
                atomList.stream().mapToDouble(atom -> atom.getCoordinates()[1]).average().getAsDouble(),
                atomList.stream().mapToDouble(atom -> atom.getCoordinates()[2]).average().getAsDouble()
        };
    }

    /**
     * Computes the root-mean-square deviation between 2 atom sets. Both atom containers must have the same number of
     * atoms.
     * @param atoms1 a collection of atoms - it is supposed to share its ordering with the 2nd container
     * @param atoms2 another collection of atoms
     * @return the RMSD value of this alignment
     * @throws IllegalArgumentException is thrown when the number of atoms does no match
     */
    public static double calculateRMSD(Stream<Atom> atoms1, Stream<Atom> atoms2) throws IllegalArgumentException {
        List<Atom> atomList1 = atoms1.collect(Collectors.toList());
        List<Atom> atomList2 = atoms2.collect(Collectors.toList());
        if(atomList1.size() != atomList2.size()) {
            throw new IllegalArgumentException("size of atom containers does not match: " + atomList1.size() + " vs " +
                    atomList2.size());
        }

        double rmsd = 0.0;
        for (int i = 0; i < atomList1.size(); i++) {
            // use squared distance for speed-up (normally the standard euclidean distance would be squared again)
            rmsd += LinearAlgebra3D.distanceFast(atomList1.get(i).getCoordinates(), atomList2.get(i).getCoordinates());
        }
        return Math.sqrt(rmsd / atomList1.size());
    }

    /**
     * Returns a stream of all residue combination whose euclidean distance is below a defined threshold.
     * @param residues a collection of residues
     * @param distanceCutoff the distance threshold below two atoms are considered to be in contact
     * @param representationScheme how to represent residues?
     * @return all residue pairs in contact according to the distance cutoff
     * @see Pair#unorderedPairsOf(List)
     * @see GroupPairDistanceCutoffFilter
     */
    public static Stream<Pair<Residue>> residuePairsInContact(Stream<Residue> residues, double distanceCutoff,
                                                              RepresentationScheme representationScheme) {
        final ResiduePairDistanceCutoffFilter distanceCutoffFilter = new ResiduePairDistanceCutoffFilter(distanceCutoff,
                representationScheme);
        return Pair.unorderedPairsOf(residues.collect(Collectors.toList())).filter(distanceCutoffFilter);
    }
}