package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

/**
 * The object representing the result of an {@link AlignmentAlgorithm}. Provides RMSD, translation vector and rotation
 * matrix necessary to recreate this particular alignment.
 * Created by S on 30.09.2016.
 */
//TODO redesign - a container also handling aligned energy profiles would be nicer
public class AlignmentResult {
    private final AtomContainer originalReference;
    private final AtomContainer originalQuery;
    private final AtomContainer alignedQuery;
    private final double alignmentScore;
    private final LinearAlgebraAtom.Transformation transformation;

    /**
     * Constructs a new alignment result container.
     * @param originalReference the original originalReference container
     * @param originalQuery the original to-be-aligned container
     * @param alignmentScore the RMSD (or any other alignment score) which was achieved for this alignment
     * @param translation the translation vector needed to recreate this alignment
     * @param rotation the rotation matrix needed to recreate this alignment
     */
    public AlignmentResult(AtomContainer originalReference, AtomContainer originalQuery, AtomContainer alignedQuery,
                           double alignmentScore, double[] translation, double[][] rotation) {
        this.originalReference = originalReference;
        this.originalQuery = originalQuery;
        this.alignedQuery = alignedQuery;
        this.alignmentScore = alignmentScore;
        this.transformation = new LinearAlgebraAtom.Transformation(translation, rotation);
    }

    /**
     * Convenience method to transform a {@link AtomContainer} based on the
     * translation and rotation described by this alignment result.
     * @param atomContainer the atom collection to transform
     * @see LinearAlgebraAtom#transform(AtomContainer, double[], double[][])
     */
    public void transform(AtomContainer atomContainer) {
        LinearAlgebraAtom.transform(atomContainer, transformation);
    }

    /**
     * Reports the RMSD achieved by this alignment.
     * @return the structural difference of both atom sets measured by the RMSD
     */
    public double getAlignmentScore() {
        return alignmentScore;
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

    public AtomContainer getOriginalReference() {
        return originalReference;
    }

    public AtomContainer getOriginalQuery() {
        return originalQuery;
    }

    public AtomContainer getAlignedQuery() {
        return alignedQuery;
    }
}
