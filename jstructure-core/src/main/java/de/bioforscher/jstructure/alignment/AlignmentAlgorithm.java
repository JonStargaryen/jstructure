package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.apache.commons.math3.linear.MatrixUtils;

/**
 * Defines the capabilities of an implementing algorithm to align two atom collections. The alignment is reported as
 * {@link StructureAlignmentResult} object which can be used to recreate the alignment - even though the coordinates of the 2nd
 * atom set handed to the {@link AlignmentAlgorithm#align(AtomContainer, AtomContainer)} method are by contract modified
 * to match the reference atom set in the best possible way.
 * TODO maybe the algorithms need some abstracted way to report whether container have to match in size - the current way does not seem to convincing
 * Created by S on 30.09.2016.
 */
public interface AlignmentAlgorithm<R extends AlignmentResult> {
    String TRANSFORMATION = "TRANSFORMATION";

    /**
     * The translation vector describing no shift whatsoever. 3D vector of zeros.
     */
    double[] NEUTRAL_TRANSLATION = new double[3];
    /**
     * The rotation matrix describing no rotation whatsoever. <tt>3 x 3</tt> identity matrix.
     */
    double[][] NEUTRAL_ROTATION = MatrixUtils.createRealIdentityMatrix(3).getData();

    /**
     * Aligns 2 atom sets.
     * @param atomContainer1 a collection of 'reference' atoms, the other collection will be superimposed onto this
     *                       select
     * @param atomContainer2 another collection of atoms whose coordinates will be manipulated, so their coordinates are
     *               optimally superimposed
     * @return statistics on the alignment
     */
    R align(AtomContainer atomContainer1, AtomContainer atomContainer2);
}
