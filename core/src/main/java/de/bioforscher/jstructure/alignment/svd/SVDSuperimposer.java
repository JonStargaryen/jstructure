package de.bioforscher.jstructure.alignment.svd;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;

/**
 * Implementation of the singular value decomposition rigid body alignment algorithm.
 * Created by S on 30.09.2016.
 */
public class SVDSuperimposer extends AbstractAlignmentAlgorithm {
    private static final Logger logger = LoggerFactory.getLogger(SVDSuperimposer.class);
    //TODO move to a more sophisticated predicate
    private final boolean backboneOnly;

    public SVDSuperimposer() {
        this(false);
    }

    public SVDSuperimposer(boolean backboneOnly) {
        this.backboneOnly = backboneOnly;
    }

    @Override
    public AlignmentResult align(AtomContainer reference, AtomContainer candidate) {
        AtomContainer originalReference = reference;
        AtomContainer originalCandidate = candidate;

        Pair<AtomContainer, AtomContainer> atomContainerPair =
                LinearAlgebraAtom.comparableAtomContainerPair(reference, candidate, backboneOnly);

        reference = atomContainerPair.getLeft();
        candidate = atomContainerPair.getRight();

        // compute centroids
        double[] centroid1 = LinearAlgebraAtom.centroid(reference);
        double[] centroid2 = LinearAlgebraAtom.centroid(candidate);
        // center atoms
        LinearAlgebraAtom.shift(reference, LinearAlgebra3D.multiply(centroid1, -1.0));
        LinearAlgebraAtom.shift(candidate, LinearAlgebra3D.multiply(centroid2, -1.0));

        // compose covariance matrix and calculate SVD
        RealMatrix matrix1 = LinearAlgebraAtom.convertToMatrix(reference);
        RealMatrix matrix2 = LinearAlgebraAtom.convertToMatrix(candidate);
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

        // compute translation
        double[] translation = LinearAlgebra3D.subtract(centroid1, LinearAlgebra3D.multiply(centroid2, rotation));

        logger.trace("rotation matrix\n{}\ntranslation vector\n{}", Arrays.deepToString(rotationMatrix.getData()),
                Arrays.toString(translation));

        /* transform 2nd atom select - employ neutral translation (3D vector of zeros), because the atoms are already
        * centered and calculate RMSD */
        LinearAlgebraAtom.transform(candidate, new LinearAlgebraAtom.Transformation(rotation));
        double rmsd = LinearAlgebraAtom.calculateRmsd(reference, candidate);

        // return alignment
        return new AlignmentResult(originalReference, originalCandidate, candidate, rmsd, translation, rotation);
    }
}
