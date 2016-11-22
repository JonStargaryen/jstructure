package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.CoordinateUtils;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

/**
 * The object representing the result of an {@link AlignmentAlgorithm}. Provides RMSD, translation vector and rotation
 * matrix necessary to recreate this particular alignment.
 * Created by S on 30.09.2016.
 */
public class AlignmentResult {
    private final AtomContainer reference;
    private final AtomContainer query;
    private final double rmsd;
    private final CoordinateUtils.Transformation transformation;

    /**
     * Constructs a new alignment result container.
     * @param reference the original reference container
     * @param query the original to-be-aligned container
     * @param rmsd the RMSD which was achieved for this alignment
     * @param translation the translation vector needed to recreate this alignment
     * @param rotation the rotation matrix needed to recreate this alignment
     */
    public AlignmentResult(AtomContainer reference, AtomContainer query, double rmsd,
                           double[] translation, double[][] rotation) {
        this.reference = reference;
        this.query = query;
        this.rmsd = rmsd;
        this.transformation = new CoordinateUtils.Transformation(translation, rotation);
    }

    /**
     * Convenience method to transform a {@link AtomContainer} based on the
     * translation and rotation described by this alignment result.
     * @param atomContainer the atom collection to transform
     * @see CoordinateUtils#transform(AtomContainer, double[], double[][])
     */
    public void transform(AtomContainer atomContainer) {
        CoordinateUtils.transform(atomContainer, transformation);
    }

    /**
     * Reports the RMSD achieved by this alignment.
     * @return the structural difference of both atom sets measured by the RMSD
     */
    public double getRmsd() {
        return rmsd;
    }

    /**
     * Access to the translation vector which can reproduce this alignment.
     * @return a 3D vector
     */
    public double[] getTranslation() {
        return transformation.getTranslation();
    }

    /**
     * Access to the rotation vector which can reproduce this alignment.
     * @return a <tt>3 x 3</tt> rotation matrix
     */
    public double[][] getRotation() {
        return transformation.getRotation();
    }

    public AtomContainer getReference() {
        return reference;
    }

    public AtomContainer getQuery() {
        return query;
    }
}
