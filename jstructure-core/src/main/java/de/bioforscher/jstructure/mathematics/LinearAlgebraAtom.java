package de.bioforscher.jstructure.mathematics;

import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A collection of commonly used function to manipulate atom coordinates. In contrast to the the functions in
 * {@link LinearAlgebra3D} most of the function provided here are somewhat high-level and are not designed to deal with
 * <code>double[]</code> but rather with {@link Atom} objects. Also they will directly modify {@link Atom} instances
 * handed over to these functions.
 * @see LinearAlgebra3D
 * Created by S on 30.09.2016.
 */
public class LinearAlgebraAtom {
    static final Logger logger = LoggerFactory.getLogger(LinearAlgebraAtom.class);

    private LinearAlgebraAtom() {
        // deny instantiation
    }

    /**
     * Computes the center of mass of a collection of atoms. In contrast to the centroid, the center of mass considers
     * the unique mass of each atom in the structure (i.e. a {@link de.bioforscher.jstructure.model.structure.Element#H})
     * has a different mass than say {@link de.bioforscher.jstructure.model.structure.Element#C} and, thus, will affect
     * the center of mass computation in a less impacting way.
     * @param atomContainer a collection of atoms
     * @return the center of mass
     * @see LinearAlgebraAtom#centroid(AtomContainer)
     */
    public static double[] centerOfMass(AtomContainer atomContainer) {
        return atomContainer.atoms().collect(StructureCollectors.toCenterOfMass());
    }

    public static TorsionAngles computeTorsionAngles(Group group1, Group group2) {
        // ensure both groups are consecutive
        if(!group1.getParentChain().equals(group2.getParentChain()) ||
                group1.getResidueNumber() + 1 != group2.getResidueNumber()) {
            throw new IllegalArgumentException("2 groups must be connected to compute torsion angles" +
                    System.lineSeparator() + "found: " + group1.getIdentifier() + " and " + group2.getIdentifier());
        }

        Atom n1 = group1.select()
                .backboneNitrogenAtoms()
                .asAtom();
        Atom ca1 = group1.select()
                .alphaCarbonAtoms()
                .asAtom();
        Atom c1 = group1.select()
                .backboneCarbonAtoms()
                .asAtom();

        Atom n2 = group2.select()
                .backboneNitrogenAtoms()
                .asAtom();
        Atom ca2 = group2.select()
                .alphaCarbonAtoms()
                .asAtom();
        Atom c2 = group2.select()
                .backboneCarbonAtoms()
                .asAtom();

        double phi = LinearAlgebra3D.torsionAngle(c1.getCoordinates(), n2.getCoordinates(), ca2.getCoordinates(), c2.getCoordinates());
        double psi = LinearAlgebra3D.torsionAngle(n1.getCoordinates(), ca1.getCoordinates(), c1.getCoordinates(), n2.getCoordinates());
        double omega = LinearAlgebra3D.torsionAngle(ca1.getCoordinates(), c1.getCoordinates(), n2.getCoordinates(), ca2.getCoordinates());

        return new TorsionAngles(phi, psi, omega);
    }

    public static class TorsionAngles {
        double phi, psi, omega;

        TorsionAngles(double phi, double psi, double omega) {
            this.phi = phi;
            this.psi = psi;
            this.omega = omega;
        }

        public double getPhi() {
            return phi;
        }

        public double getPsi() {
            return psi;
        }

        public double getOmega() {
            return omega;
        }
    }

    /**
     * Computes the maximal extent of this protein in any given spatial direction to the centroid of this structure.
     * @param atomContainer a collection of atoms
     * @return the maximal distance occurring between the centroid and any other atom
     * @see LinearAlgebraAtom#maximalExtent(AtomContainer, double[])
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
        return atomContainer.select()
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

    public static class Transformation extends FeatureContainerEntry {
        private double[] translation;
        private double[][] rotation;

        /**
         * Mere shift.
         * @param translation the translation vector
         */
        public Transformation(double[] translation) {
            super(null);
            this.translation = translation;
        }

        /**
         * Rototranslation.
         * @param translation the translation vector
         * @param rotation the rotation matrix
         */
        public Transformation(double[] translation, double[][] rotation) {
            super(null);
            this.translation = translation;
            this.rotation = rotation;
        }

        /**
         * Mere rotation.
         * @param rotation the rotation matrix
         */
        public Transformation(double[][] rotation) {
            super(null);
            this.rotation = rotation;
        }

        public Atom transformCoordinates(Atom atom) {
            double[] vector = atom.getCoordinates();
            logger.trace("initial atom {}", atom.getPdbRepresentation());

            // apply rotation if needed
            if(rotation != null) {
                atom.setCoordinates(LinearAlgebra3D.multiply(vector, rotation));
            }

            vector = atom.getCoordinates();
            // apply transformation if needed
            if(translation != null) {
                atom.setCoordinates(LinearAlgebra3D.add(vector, translation));
            }

            logger.trace("transf. atom {}", atom.getPdbRepresentation());
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
     * @see LinearAlgebraAtom#centroid(AtomContainer)
     * @see LinearAlgebraAtom#shift(AtomContainer, double[])
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
        Pair<GroupContainer, GroupContainer> groupContainerPair = comparableGroupContainerPair(atomContainer1, atomContainer2);
        double msd = Combinatorics.sequentialPairsOf(groupContainerPair.getLeft().getAtoms(), groupContainerPair.getRight().getAtoms())
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
     * @param container1 a collection of reference atoms
     * @param container2 a collection of candidate atoms
     * @param minimalSetOfAtomNames the lower bound of required atom names present (e.g. by setting CA, every time an
     *                              amino acid missing the alpha carbon will result in an exception thrown)
     * @param maximalSetOfAtomNames the upper bound of required atom names present (e.g. by setting CA, everything else
     *                              will be dropped, even when both amino acids would share more atoms)
     * @return a pair of both collections which can now be aligned
     */
    public static Pair<GroupContainer, GroupContainer> comparableGroupContainerPair(AtomContainer container1,
                                                                                    AtomContainer container2,
                                                                                    Set<String> minimalSetOfAtomNames,
                                                                                    Set<String> maximalSetOfAtomNames) {
        GroupContainer groupContainer1 = cloneIntoGroupContainer(container1);
        GroupContainer groupContainer2 = cloneIntoGroupContainer(container2);

        int limitingSize = Math.min(groupContainer1.getGroups().size(), groupContainer2.getGroups().size());

        List<Group> groups1 = new ArrayList<>();
        List<Group> groups2 = new ArrayList<>();
        for(int groupIndex = 0; groupIndex < limitingSize; groupIndex++) {
            Group group1 = groupContainer1.getGroups().get(groupIndex);
            Group group2 = groupContainer2.getGroups().get(groupIndex);
            Pair<List<Atom>, List<Atom>> sharedAtoms = selectSharedAtoms(group1,
                    group2,
                    minimalSetOfAtomNames,
                    maximalSetOfAtomNames);

            // remove additional/not-shared atoms
            group1.getAtoms().clear();
            group1.getAtoms().addAll(sharedAtoms.getLeft());
            group2.getAtoms().clear();
            group2.getAtoms().addAll(sharedAtoms.getRight());

            logger.trace("shared atoms between {} and {}: {}", group1, group2, sharedAtoms);

            groups1.add(group1);
            groups2.add(group2);
        }

        return new Pair<>(new Chain(groups1), new Chain(groups2));
    }

    /**
     * Returns the set of atoms shared by both containers.
     * @see LinearAlgebraAtom#comparableGroupContainerPair(AtomContainer, AtomContainer, Set, Set)
     */
    public static Pair<GroupContainer, GroupContainer> comparableGroupContainerPair(AtomContainer atomContainer1,
                                                                                    AtomContainer atomContainer2) {
        return comparableGroupContainerPair(atomContainer1,
                atomContainer2,
                Collections.emptySet(),
                Collections.emptySet());
    }

    private static Pair<List<Atom>, List<Atom>> selectSharedAtoms(Group group1,
                                                                  Group group2,
                                                                  Set<String> minimalSetOfAtomNames,
                                                                  Set<String> maximalSetOfAtomNames) {
        // no ordering is enforced - rather the filtering by atom names ensures the identical occurrence of atoms
        Set<String> sharedAtomNames = determineSharedAtomNames(group1,
                group2,
                minimalSetOfAtomNames,
                maximalSetOfAtomNames);

        return new Pair<>(selectAtoms(group1, sharedAtomNames), selectAtoms(group2, sharedAtomNames));
    }

    private static List<Atom> selectAtoms(Group group, Set<String> atomNamesToSelect) {
        // fix 02/08/17 - mere retaining does not ensure correct ordering
        return atomNamesToSelect.stream()
                .flatMap(name -> group.atoms()
                        .filter(atom -> atom.getName().equals(name)))
                .collect(Collectors.toList());
    }

    private static Set<String> determineSharedAtomNames(Group group1,
                                                        Group group2,
                                                        Set<String> minimalSetOfAtomNames,
                                                        Set<String> maximalSetOfAtomNames) {
        // present atom names in each group, determine the shared ones
        Set<String> sharedAtomNames = group1.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet());
        sharedAtomNames.retainAll(group2.atoms()
                .map(Atom::getName)
                .collect(Collectors.toSet()));

        // fail if the minimal set of atoms is not fulfilled
        if(!sharedAtomNames.containsAll(minimalSetOfAtomNames)) {
            throw new IllegalArgumentException("alignment could not fulfill minimal required atom names" +
                    System.lineSeparator() + group1.getIdentifier() + " atoms: " + group1.atoms().map(Atom::getName).collect(Collectors.joining(", ")) +
                    System.lineSeparator() + group2.getIdentifier() + " atoms: " + group2.atoms().map(Atom::getName).collect(Collectors.joining(", ")) +
                    System.lineSeparator() + "shared: " + sharedAtomNames.stream().collect(Collectors.joining(", ")) +
                    System.lineSeparator() + "required: " + minimalSetOfAtomNames.stream().collect(Collectors.joining(", ")));
        }

        // ignore all atoms not explicitly wanted in the alignment
        if(!maximalSetOfAtomNames.isEmpty()) {
            sharedAtomNames.retainAll(maximalSetOfAtomNames);
        }

        return sharedAtomNames;
    }

    /**
     * Maps an {@link AtomContainer} to a distinct copy of itself.
     * @param atomContainer the container to clone
     * @return a {@link GroupContainer} (i.e. motivated by interfaces handling {@link AtomContainer}) which can readily
     * manipulated without affecting the original entry
     */
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