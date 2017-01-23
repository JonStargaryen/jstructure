package de.bioforscher.jstructure.alignment;

import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;

import java.util.Collections;
import java.util.Set;

/**
 * The abstract implementation of alignment algorithms.
 * Created by S on 10.11.2016.
 */
public abstract class AbstractAlignmentAlgorithm implements AlignmentAlgorithm {
    final Set<String> minimalSetOfAtomNames;
    final Set<String> maximalSetOfAtomNames;
    public static final String FRAGMENT_RMSD = "FRAGMENT_RMSD";

    public AbstractAlignmentAlgorithm() {
        this(Collections.emptySet(), AminoAcidFamily.ATOM_NAMES.ALL_ATOM_NAMES);
    }

    public AbstractAlignmentAlgorithm(Set<String> minimalSetOfAtomNames, Set<String> maximalSetOfAtomNames) {
        this.minimalSetOfAtomNames = minimalSetOfAtomNames;
        this.maximalSetOfAtomNames = maximalSetOfAtomNames;
    }
}
