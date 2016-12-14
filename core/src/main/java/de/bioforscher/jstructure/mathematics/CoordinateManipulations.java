package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import static de.bioforscher.jstructure.model.structure.selection.Selection.on;

/**
 * A collection of commonly used function to manipulate atom coordinates. In contrast to the the functions in
 * {@link LinearAlgebra3D} most of the function provided here are somewhat high-level and are not designed to deal with
 * <code>double[]</code> but rather with {@link Atom} objects. Also they will directly modify {@link Atom} instances
 * handed over to these functions.
 * @see LinearAlgebra3D
 * Created by S on 30.09.2016.
 */
public class CoordinateManipulations {
    static final Logger logger = LoggerFactory.getLogger(CoordinateManipulations.class);

    private CoordinateManipulations() {
        // deny instantiation
    }

    /**
     * Computes the center of mass of a collection of atoms. In contrast to the centroid, the center of mass considers
     * the unique mass of each atom in the structure (i.e. a {@link de.bioforscher.jstructure.model.structure.Element#H})
     * has a different mass than say {@link de.bioforscher.jstructure.model.structure.Element#C} and, thus, will affect
     * the center of mass computation in a less impactful way.
     * @param atomContainer a collection of atoms
     * @return the center of mass
     * @see CoordinateManipulations#centroid(AtomContainer)
     */
    public static double[] centerOfMass(AtomContainer atomContainer) {
        return atomContainer.atoms().collect(StructureCollectors.toCenterOfMass());
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atomContainer a collection of atoms
     * @return the maximal distance occurring between the centroid and any other atom
     * @see CoordinateManipulations#maximalExtent(AtomContainer, double[])
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
        return on(atomContainer)
                .alphaCarbonAtoms()
                .asFilteredAtoms()
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
        transform(atomContainer, new Transformation(translation, rotation));
    }

    /**
     * Employs a transformation (i.e. a rotation & translation or 'rototranslation') on a collection of atoms. The
     * operation will manipulate the provided atom's coordinates directly (rather than creating new Objects).
     * @param atomContainer a collection of atoms
     * @param transformation the transformation to employ
     */
    public static void transform(AtomContainer atomContainer, Transformation transformation) {
        atomContainer.atoms().forEach(transformation::transformCoordinates);
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
            logger.trace("initial atom {}", atom.composePDBRecord());

            // apply rotation if needed
            if(rotation != null) {
                atom.setCoordinates(LinearAlgebra3D.multiply(vector, rotation));
            }

            vector = atom.getCoordinates();
            // apply transformation if needed
            if(translation != null) {
                atom.setCoordinates(LinearAlgebra3D.add(vector, translation));
            }

            logger.trace("transf. atom {}", atom.composePDBRecord());
            return atom;
        }

        public double[] getTranslation() {
            return translation;
        }

        public double[][] getRotation() {
            return rotation;
        }
    }

    /**
     * Centers a collection of atoms. This is achieved by computing the centroid of these points and subtracting this
     * point from each atom's coordinates. The operation will manipulate the provided atom's coordinates directly
     * (rather than creating new Objects).
     * @param atomContainer the collection of atoms to center
     * @see CoordinateManipulations#centroid(AtomContainer)
     * @see CoordinateManipulations#shift(AtomContainer, double[])
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
        return atomContainer.atoms().collect(StructureCollectors.toCentroid());
    }

    /**
     * Computes the root-mean-square deviation between 2 atom sets.
     * @param atomContainer1 a collection of atoms
     * @param atomContainer2 another collection of atoms
     * @return the RMSD value of this alignment
     */
    public static double calculateRmsd(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        Pair<AtomContainer, AtomContainer> atomContainerPair = comparableAtomContainerPair(atomContainer1, atomContainer2);
        double msd = Combinatorics.sequentialPairsOf(atomContainerPair.getLeft().getAtoms(), atomContainerPair.getRight().getAtoms())
                .mapToDouble(pair -> LinearAlgebra3D.distanceFast(pair.getLeft().getCoordinates(), pair.getRight().getCoordinates()))
                .average()
                .orElseThrow(() -> new IllegalArgumentException("cannot compute rmsd for empty or non-intersecting containers"));
        return Math.sqrt(msd);
    }

    /**
     * Creates a comparable entity of 2 {@link AtomContainer} objects. They fulfill certain criteria:
     * <ul>
     *     <li>{@link Atom}s are cloned
     *     <li>for each {@link Group} only shared atoms are retained</li>
     *     <li>{@link Atom}s are in identical ordering</li>
     *     <li>{@link Atom}s are wrapped in {@link Group} objects if they were initially</li>
     * </ul>
     * @param atomContainer1 a collection of reference atoms
     * @param atomContainer2 a collection of candidate atoms
     * @return a pair of both collections which can now be aligned
     */
    public static Pair<AtomContainer, AtomContainer> comparableAtomContainerPair(AtomContainer atomContainer1,
                                                                                 AtomContainer atomContainer2,
                                                                                 boolean backboneOnly) {
        GroupContainer groupContainer1 = cloneIntoGroupContainer(atomContainer1);
        GroupContainer groupContainer2 = cloneIntoGroupContainer(atomContainer2);

        if(groupContainer1.getGroups().size() != groupContainer2.getGroups().size()) {
            throw new IllegalArgumentException("cannot compare atom containers of different size: " + atomContainer1.getIdentifier() + " : " + atomContainer2.getIdentifier());
        }

        List<Group> groups1 = new ArrayList<>();
        List<Group> groups2 = new ArrayList<>();
        for(int groupIndex = 0; groupIndex < groupContainer1.getGroups().size(); groupIndex++) {
            Group group1 = groupContainer1.getGroups().get(groupIndex);
            Group group2 = groupContainer2.getGroups().get(groupIndex);
            Pair<List<Atom>, List<Atom>> sharedAtoms = selectSharedAtoms(group1, group2, backboneOnly);

            // remove additional/not-shared atoms
            group1.getAtoms().retainAll(sharedAtoms.getLeft());
            group2.getAtoms().retainAll(sharedAtoms.getRight());

            logger.trace("shared atoms between {} and {}: {}", group1, group2, sharedAtoms);
            groups1.add(group1);
            groups2.add(group2);
        }

        return new Pair<>(new Chain(groups1), new Chain(groups2));
    }

    public static Pair<AtomContainer, AtomContainer> comparableAtomContainerPair(AtomContainer atomContainer1,
                                                                                 AtomContainer atomContainer2) {
        return comparableAtomContainerPair(atomContainer1, atomContainer2, false);
    }

    private static Pair<List<Atom>, List<Atom>> selectSharedAtoms(Group group1, Group group2, boolean backboneOnly) {
        // no ordering is enforced - rather the filtering by atom names ensures the identical occurrence of atoms
        Set<String> sharedAtomNames = determineSharedAtomNames(group1, group2, backboneOnly);
        return new Pair<>(group1.atoms()
                .filter(atom -> sharedAtomNames.contains(atom.getName()))
                .collect(Collectors.toList()),
                group2.atoms()
                .filter(atom -> sharedAtomNames.contains(atom.getName()))
                .collect(Collectors.toList()));
    }

    private static Set<String> determineSharedAtomNames(Group group1, Group group2, boolean backboneOnly) {
        // present atom names in each group, determine the shared ones
        Set<String> sharedAtomNames = group1.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet());
        sharedAtomNames.retainAll(group2.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet()));

        if(backboneOnly) {
            sharedAtomNames.retainAll(AminoAcid.ATOM_NAMES.BACKBONE_ATOM_NAMES);
        }

        return sharedAtomNames;
    }

    private static GroupContainer cloneIntoGroupContainer(AtomContainer atomContainer) {
        return atomContainer.atoms()
                // map to group level - will still collect 'dangling' atoms into a group
                .map(Atom::getParentGroup)
                .distinct()
                // clone
                .map(Group::new)
                .collect(StructureCollectors.toGroupContainer());
    }
}