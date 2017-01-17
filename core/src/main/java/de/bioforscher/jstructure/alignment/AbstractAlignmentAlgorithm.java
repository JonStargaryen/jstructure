package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collections;
import java.util.Set;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    final static Logger logger = LoggerFactory.getLogger(AbstractAlignmentAlgorithm.class);
    final Set<String> minimalSetOfAtomNames;
    final Set<String> maximalSetOfAtomNames;
    public static final String FRAGMENT_RMSD = "FRAGMENT_RMSD";

    AbstractAlignmentAlgorithm() {
        this(Collections.emptySet(), AminoAcidFamily.ATOM_NAMES.ALL_ATOM_NAMES);
    }

    AbstractAlignmentAlgorithm(Set<String> minimalSetOfAtomNames, Set<String> maximalSetOfAtomNames) {
        this.minimalSetOfAtomNames = minimalSetOfAtomNames;
        this.maximalSetOfAtomNames = maximalSetOfAtomNames;
    }
}
