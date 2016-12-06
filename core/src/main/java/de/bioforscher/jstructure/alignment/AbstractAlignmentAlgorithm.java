package de.bioforscher.jstructure.alignment;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    final static Logger logger = LoggerFactory.getLogger(AbstractAlignmentAlgorithm.class);

    public enum FeatureNames {
        FRAGMENT_RMSD
    }
}
