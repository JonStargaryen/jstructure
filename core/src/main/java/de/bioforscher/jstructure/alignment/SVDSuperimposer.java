package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.AtomContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Implementation of the singular value decomposition rigid body alignment algorithm.
 * Created by S on 30.09.2016.
 */
public class SVDSuperimposer implements AlignmentAlgorithm {
    public static final AtomNameFilter DEFAULT_ALIGNED_ATOMS = AtomNameFilter.CA_ATOM_FILTER;
    private AtomNameFilter alignedAtomNameFilter;

    public SVDSuperimposer(AtomNameFilter alignedAtomNameFilter) {
        this.alignedAtomNameFilter = alignedAtomNameFilter;
    }

    /**
     * The default constructor will filter for alpha carbons and align them.
     */
    public SVDSuperimposer() {
        this(DEFAULT_ALIGNED_ATOMS);
    }

    @Override
    public AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        List<Atom> atomList1 = atomContainer1.atoms().filter(alignedAtomNameFilter).collect(Collectors.toList());
        List<Atom> atomList2 = atomContainer2.atoms().filter(alignedAtomNameFilter).collect(Collectors.toList());

        // compute centroids
        double[] centroid1 = LinearAlgebra3D.multiply(CoordinateUtils.centroid(atomList1.stream()), -1.0);
        double[] centroid2 = LinearAlgebra3D.multiply(CoordinateUtils.centroid(atomList2.stream()), -1.0);
        // center atoms
        CoordinateUtils.shift(atomList1.stream(), centroid1);
        CoordinateUtils.shift(atomList2.stream(), centroid2);

        // compose covariance matrix and calculate SVD
        RealMatrix matrix1 = CoordinateUtils.convertToMatrix(atomList1.stream());
        RealMatrix matrix2 = CoordinateUtils.convertToMatrix(atomList2.stream());
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
        RealMatrix referenceCentroidMatrix = new Array2DRowRealMatrix(centroid1).transpose();
        RealMatrix fragmentCentroidMatrix = new Array2DRowRealMatrix(centroid2).transpose();
        double[] translation = referenceCentroidMatrix.subtract(
                fragmentCentroidMatrix.multiply(rotationMatrix)).getRow(0);

        /* transform 2nd atom stream - employ neutral translation (3D vector of zeros), because the atoms are already
        * centered and calculate RMSD */
        double rmsd = CoordinateUtils.calculateRMSD(atomList1.stream(),
                CoordinateUtils.transform(atomList2.stream(),
                NEUTRAL_TRANSLATION,
                rotation));
        // return alignment
        return new AlignmentResult(rmsd, translation, rotation);
    }

    @Override
    public AtomNameFilter getAlignedAtomNameFilter() {
        return alignedAtomNameFilter;
    }
}
