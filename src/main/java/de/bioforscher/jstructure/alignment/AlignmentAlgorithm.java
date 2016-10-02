package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.AtomContainer;
import org.apache.commons.math3.linear.MatrixUtils;

/**
 * Defines the capabilities of an implementing algorithm to align two atom collections. The alignment is reported as
 * {@link AlignmentResult} object which can be used to recreate the alignment - even though the coordinates of the 2nd
 * atom set handed to the {@link AlignmentAlgorithm#align(AtomContainer, AtomContainer)} method are by contract modified
 * to match the reference atom set in the best possible way.
 * TODO maybe the algorithms need some abstracted way to report what atoms are aligned and whether container have to match in size - the current way does not seem to convincing
 * Created by S on 30.09.2016.
 */
public interface AlignmentAlgorithm {
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
     * @param atomContainer1 a collection of 'reference' atoms, the other collection will be superimposed on this stream
     * @param atomContainer2 another collection of atoms whose coordinates will be manipulated, so their coordinates are
     *               optimally superimposed
     * @return statistics on the alignment
     */
    AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2);

    /**
     * Reports what atoms are aligned by this algorithm. Some will align all atoms handed to them, some will select the
     * considered atoms by name.
     * @return the {@link AtomContainer.AtomNameFilter} the implementing
     * algorithm applies
     */
    AtomContainer.AtomNameFilter getAlignedAtomNameFilter();
}
