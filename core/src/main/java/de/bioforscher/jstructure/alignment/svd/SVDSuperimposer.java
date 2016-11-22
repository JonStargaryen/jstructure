package de.bioforscher.jstructure.alignment.svd;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Implementation of the singular value decomposition rigid body alignment algorithm.
 * Created by S on 30.09.2016.
 */
public class SVDSuperimposer extends AbstractAlignmentAlgorithm {
    final Logger logger = LoggerFactory.getLogger(SVDSuperimposer.class);
    public static final List<String> DEFAULT_ALIGNED_ATOMS = Collections.singletonList(AminoAcid.ATOM_NAMES.CA_ATOM_NAME);

    public SVDSuperimposer(List<String> alignedAtomNames) {
        super(alignedAtomNames);
    }

    /**
     * The default constructor will select for alpha carbons and align them.
     */
    public SVDSuperimposer() {
        this(DEFAULT_ALIGNED_ATOMS);
    }

    @Override
    protected AlignmentResult alignInternal(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        //TODO move to abstract impl and provide flag/interface
        if(atomContainer1.getAtoms().size() != atomContainer2.getAtoms().size()) {
            logger.error("arrays do not match in size\n{}\n{}",
                    atomContainer1.composePDBRecord(),
                    atomContainer2.composePDBRecord());
            throw new DimensionMismatchException(atomContainer1.getAtoms().size(), atomContainer2.getAtoms().size());
        }

        // compute centroids
        double[] centroid1 = CoordinateUtils.centroid(atomContainer1);
        double[] centroid2 = CoordinateUtils.centroid(atomContainer2);
        // center atoms
        CoordinateUtils.shift(atomContainer1, LinearAlgebra3D.multiply(centroid1, -1.0));
        CoordinateUtils.shift(atomContainer2, LinearAlgebra3D.multiply(centroid2, -1.0));

        // compose covariance matrix and calculate SVD
        RealMatrix matrix1 = CoordinateUtils.convertToMatrix(atomContainer1);
        RealMatrix matrix2 = CoordinateUtils.convertToMatrix(atomContainer2);
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
        CoordinateUtils.transform(atomContainer2, new CoordinateUtils.Transformation(rotation));
        double rmsd = CoordinateUtils.calculateRMSD(atomContainer1, atomContainer2);
        // return alignment
        return new AlignmentResult(atomContainer1, atomContainer2, rmsd, translation, rotation);
    }
}
