package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.List;

/**
 * The new take on structural alignments.
 * Created by bittrich on 6/19/17.
 */
public class StructureAligner {
    private static final Logger logger = LoggerFactory.getLogger(StructureAligner.class);
    // original containers
    private final GroupContainer referenceOriginal;
    private final GroupContainer queryOriginal;
    // potentially cloned (i.e. decoupled containers)
    private final GroupContainer reference;
    private final GroupContainer query;
    // actually selected atoms used for alignment of both containers
    private final AtomContainer referenceSelectedAtoms;
    private final AtomContainer querySelectedAtoms;
    private final List<Pair<Atom, Atom>> atomMapping;
    private final AlignmentPolicy.ManipulationBehavior manipulationBehavior;

    StructureAligner(StructureAlignerBuilder builder) {
        this.referenceOriginal = builder.reference;
        this.queryOriginal = builder.query;

        this.reference = referenceOriginal.createCopy();
        this.query = queryOriginal.createCopy();

        // determine mapping
        this.atomMapping = builder.atomMapping.determineAtomMapping(reference, query);
        this.referenceSelectedAtoms = atomMapping.stream()
                .map(Pair::getLeft)
                .collect(StructureCollectors.toAtomContainer());
        this.querySelectedAtoms = atomMapping.stream()
                .map(Pair::getRight)
                .collect(StructureCollectors.toAtomContainer());

        this.manipulationBehavior = builder.manipulationBehavior;
    }

    public Alignment align() {
        // calculate centroids and center atoms
        double[] centroid1 = referenceSelectedAtoms.calculate().center().getValue();
        double[] centroid2 = querySelectedAtoms.calculate().center().getValue();

        // compose covariance matrix and calculate SVD
        RealMatrix matrix1 = convertToMatrix(referenceSelectedAtoms);
        RealMatrix matrix2 = convertToMatrix(querySelectedAtoms);
        RealMatrix covariance = matrix2.transpose().multiply(matrix1);
        SingularValueDecomposition svd = new SingularValueDecomposition(covariance);
        // R = (V * U')'
        RealMatrix ut = svd.getU().transpose();
        RealMatrix rotationMatrix = svd.getV().multiply(ut).transpose();
        // check if reflection
        if(new LUDecomposition(rotationMatrix).getDeterminant() < 0) {
            RealMatrix v = svd.getV().transpose();
            v.setEntry(2, 0, (0 - v.getEntry(2, 0)));
            v.setEntry(2, 1, (0 - v.getEntry(2, 1)));
            v.setEntry(2, 2, (0 - v.getEntry(2, 2)));
            rotationMatrix = v.transpose().multiply(ut).transpose();
        }
        double[][] rotation = rotationMatrix.getData();

        // calculate translation
        double[] translation = LinearAlgebra.on(centroid1).subtract(LinearAlgebra.on(centroid2).multiply(rotation)).getValue();

        logger.trace("rotation matrix\n{}\ntranslation vector\n{}", Arrays.deepToString(rotationMatrix.getData()),
                Arrays.toString(translation));

        /* transform 2nd atom select - employ neutral translation (3D vector of zeros), because the atoms are already
        * centered and calculate RMSD */
        querySelectedAtoms.calculate().transform(new Transformation(rotation));
        double rmsd = calculateRmsd();

        Transformation transformation = new Transformation(translation, rotation);
        // superimpose query onto reference
        query.calculate().transform(transformation);
        if(manipulationBehavior == AlignmentPolicy.ManipulationBehavior.INPLACE) {
            // if requested: do the same for the original query
            queryOriginal.calculate().transform(transformation);
        }

        // return alignment
        return new Alignment(referenceOriginal,
                queryOriginal,
                query,
                transformation,
                rmsd);
    }

    /**
     * Computes the root-mean square deviation between 2 sets of atoms.
     * @return the RMSD value of the alignment
     * @throws IllegalArgumentException if no matching atom pairs were provided
     */
    private double calculateRmsd() {
        double msd = atomMapping.stream()
                .mapToDouble(pair -> LinearAlgebra.on(pair.getLeft().getCoordinates()).distanceFast(pair.getRight().getCoordinates()))
                .average()
                .orElseThrow(() -> new IllegalArgumentException("cannot calculate rmsd for empty or non-intersecting containers"));
        return Math.sqrt(msd);
    }

    /**
     * Converts a collection of atoms to a <tt>n x 3</tt> matrix, where <tt>n</tt> is equal to the number of processed
     * atoms.
     * @param atomContainer a collection of atoms
     * @return a matrix containing the coordinates of all atoms
     */
    private RealMatrix convertToMatrix(AtomContainer atomContainer) {
        List<Atom> atomList = atomContainer.getAtoms();
        double[][] matrix = new double[atomList.size()][3];
        for (int index = 0; index < atomList.size(); index++) {
            matrix[index] = atomList.get(index).getCoordinates();
        }
        return new Array2DRowRealMatrix(matrix);
    }

    public static MatchingBehaviorStep builder(GroupContainer reference, GroupContainer query) {
        return new MatchingBehaviorStep(reference, query);
    }

    public static class MatchingBehaviorStep {
        final GroupContainer reference;
        final GroupContainer query;
        AlignmentPolicy.AtomMapping atomMapping;

        MatchingBehaviorStep(GroupContainer reference, GroupContainer query) {
            this.reference = reference;
            this.query = query;
        }

        public ManipulationBehaviorStep matchingBehavior(AlignmentPolicy.AtomMapping atomMapping) {
            this.atomMapping = atomMapping;
            return new ManipulationBehaviorStep(this);
        }
    }

    public static class ManipulationBehaviorStep {
        final GroupContainer reference;
        final GroupContainer query;
        final AlignmentPolicy.AtomMapping atomMapping;
        AlignmentPolicy.ManipulationBehavior manipulationBehavior;

        ManipulationBehaviorStep(MatchingBehaviorStep matchingBehaviorStep) {
            this.reference = matchingBehaviorStep.reference;
            this.query = matchingBehaviorStep.query;
            this.atomMapping = matchingBehaviorStep.atomMapping;
        }

        public StructureAlignerBuilder manipulationBehavior(AlignmentPolicy.ManipulationBehavior manipulationBehavior) {
            this.manipulationBehavior = manipulationBehavior;
            return new StructureAlignerBuilder(this);
        }
    }

    public static class StructureAlignerBuilder {
        final GroupContainer reference;
        final GroupContainer query;
        final AlignmentPolicy.AtomMapping atomMapping;
        final AlignmentPolicy.ManipulationBehavior manipulationBehavior;

        StructureAlignerBuilder(ManipulationBehaviorStep manipulationBehaviorStep) {
            this.reference = manipulationBehaviorStep.reference;
            this.query = manipulationBehaviorStep.query;
            this.atomMapping = manipulationBehaviorStep.atomMapping;
            this.manipulationBehavior = manipulationBehaviorStep.manipulationBehavior;
        }

        public Alignment align() {
            return new StructureAligner(this).align();
        }
    }
}
