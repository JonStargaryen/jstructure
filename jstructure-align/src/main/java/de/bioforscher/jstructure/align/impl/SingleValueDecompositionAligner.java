package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.*;
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
import java.util.stream.IntStream;

/**
 * A structure alignment algorithm based on SVD.
 * Created by S on 12.07.2017.
 */
public class SingleValueDecompositionAligner implements StructureAligner {
    private static final Logger logger = LoggerFactory.getLogger(SingleValueDecompositionAligner.class);

    @Override
    public StructureAlignmentResult align(StructureAlignmentQuery structureAlignmentQuery) throws AlignmentException {
        try {
            GroupContainer referenceOriginal = structureAlignmentQuery.getReference();
            GroupContainer queryOriginal = structureAlignmentQuery.getQuery();

            // determine mapping
            List<Pair<Atom, Atom>> atomMapping = structureAlignmentQuery.getAtomMapping()
                    .determineAtomMapping(referenceOriginal, queryOriginal);
            AtomContainer referenceSelectedAtoms = atomMapping.stream()
                    .map(Pair::getLeft)
                    .collect(StructureCollectors.toIsolatedStructure());
            AtomContainer querySelectedAtoms = atomMapping.stream()
                    .map(Pair::getRight)
                    .collect(StructureCollectors.toIsolatedStructure());

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
            if (new LUDecomposition(rotationMatrix).getDeterminant() < 0) {
                RealMatrix v = svd.getV().transpose();
                v.setEntry(2, 0, (0 - v.getEntry(2, 0)));
                v.setEntry(2, 1, (0 - v.getEntry(2, 1)));
                v.setEntry(2, 2, (0 - v.getEntry(2, 2)));
                rotationMatrix = v.transpose().multiply(ut).transpose();
            }
            double[][] rotation = rotationMatrix.getData();

            // calculate translation
            double[] translation = LinearAlgebra.on(centroid1)
                    .subtract(LinearAlgebra.on(centroid2)
                            .multiply(rotation))
                    .getValue();

            logger.trace("rotation matrix\n{}\ntranslation vector\n{}", Arrays.deepToString(rotationMatrix.getData()),
                    Arrays.toString(translation));

            // transform 2nd atom select - employ neutral translation (3D vector of zeros), because the atoms are
            // already centered and calculate RMSD
            querySelectedAtoms.calculate().transform(new Transformation(rotation));
            double rmsd = calculateRmsd(referenceSelectedAtoms, querySelectedAtoms);

            Transformation transformation = new Transformation(translation, rotation);
            // superimpose query onto reference
            querySelectedAtoms.calculate().transform(transformation);

            // return alignment
            return new StructureAlignmentResultImpl(referenceOriginal,
                    queryOriginal,
                    referenceSelectedAtoms,
                    querySelectedAtoms,
                    atomMapping,
                    transformation,
                    rmsd);
        } catch (Exception e) {
            throw new AlignmentException(e);
        }
    }

    /**
     * Computes the root-mean square deviation between 2 sets of atoms.
     * @param reference container 1 - must arrange atoms in the exact same manner
     * @param query container 2 - must arrange atoms in the exact same manner
     * @return the RMSD value of the alignment
     * @throws IllegalArgumentException if no matching atom pairs were provided
     */
    private double calculateRmsd(AtomContainer reference, AtomContainer query) {
        double msd = IntStream.range(0, reference.getAtoms().size())
                .mapToDouble(atomIndex -> LinearAlgebra.on(reference.getAtoms().get(atomIndex))
                        .distanceFast(query.getAtoms().get(atomIndex)))
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
}