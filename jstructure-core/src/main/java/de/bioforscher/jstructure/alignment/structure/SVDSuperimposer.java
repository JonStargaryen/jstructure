package de.bioforscher.jstructure.alignment.structure;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Implementation of the singular value decomposition rigid body alignment algorithm.
 * Created by S on 30.09.2016.
 */
public class SVDSuperimposer extends AbstractAlignmentAlgorithm<StructureAlignmentResult> {
    private static final Logger logger = LoggerFactory.getLogger(SVDSuperimposer.class);
    private static final Set<String> ALPHA_CARBON_NAMES = Stream.of(AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME).collect(Collectors.toSet());
    public static final SVDSuperimposer ALPHA_CARBON_SVD_INSTANCE = new SVDSuperimposer(ALPHA_CARBON_NAMES, ALPHA_CARBON_NAMES);
    public static final SVDSuperimposer BACKBONE_SVD_INSTANCE = new SVDSuperimposer(AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES,
            AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES);

    public SVDSuperimposer() {
        super();
    }

    public SVDSuperimposer(Set<String> minimalSetOfAtomNames, Set<String> maximalSetOfAtomNames) {
        super(minimalSetOfAtomNames, maximalSetOfAtomNames);
    }

    @Override
    public StructureAlignmentResult align(AtomContainer reference, AtomContainer query) {
        AtomContainer originalReference = reference;
        AtomContainer originalCandidate = query;

        Pair<GroupContainer, GroupContainer> atomContainerPair =
                LinearAlgebraAtom.comparableGroupContainerPair(reference,
                        query,
                        minimalSetOfAtomNames,
                        maximalSetOfAtomNames);

        reference = atomContainerPair.getLeft();
        query = atomContainerPair.getRight();

        // compute centroids
        double[] centroid1 = LinearAlgebraAtom.centroid(reference);
        double[] centroid2 = LinearAlgebraAtom.centroid(query);
        // center atoms
        LinearAlgebraAtom.shift(reference, LinearAlgebra3D.multiply(centroid1, -1.0));
        LinearAlgebraAtom.shift(query, LinearAlgebra3D.multiply(centroid2, -1.0));

        // compose covariance matrix and calculate SVD
        RealMatrix matrix1 = LinearAlgebraAtom.convertToMatrix(reference);
        RealMatrix matrix2 = LinearAlgebraAtom.convertToMatrix(query);
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
        LinearAlgebraAtom.transform(query, new LinearAlgebraAtom.Transformation(rotation));
        double rmsd = LinearAlgebraAtom.calculateRmsd(reference, query);

        // return alignment
        return new StructureAlignmentResult(originalReference, originalCandidate, query, rmsd, translation, rotation);
    }
}
