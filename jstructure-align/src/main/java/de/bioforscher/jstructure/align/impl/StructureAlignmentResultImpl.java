package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.StructureAlignmentResult;
import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.List;

/**
 * The result of an structure alignment.
 * Created by bittrich on 6/19/17.
 */
public class StructureAlignmentResultImpl implements StructureAlignmentResult {
    private final AtomContainer originalReference;
    private final AtomContainer originalQuery;
    private final AtomContainer alignedReference;
    private final AtomContainer alignedQuery;
    private final List<Pair<Atom, Atom>> atomPairing;
    private final Transformation transformation;
    private final double alignmentScore;

    StructureAlignmentResultImpl(AtomContainer originalReference,
                                 AtomContainer originalQuery,
                                 AtomContainer alignedReference,
                                 AtomContainer alignedQuery,
                                 List<Pair<Atom, Atom>> atomPairing,
                                 Transformation transformation,
                                 double alignmentScore) {
        this.originalReference = originalReference;
        this.originalQuery = originalQuery;
        this.alignedReference = alignedReference;
        this.alignedQuery = alignedQuery;
        this.atomPairing = atomPairing;
        this.transformation = transformation;
        this.alignmentScore = alignmentScore;
    }

    @Override
    public AtomContainer getOriginalReference() {
        return originalReference;
    }

    @Override
    public AtomContainer getOriginalQuery() {
        return originalQuery;
    }

    @Override
    public AtomContainer getAlignedReference() {
        return alignedReference;
    }

    @Override
    public AtomContainer getAlignedQuery() {
        return alignedQuery;
    }

    @Override
    public List<Pair<Atom, Atom>> getAtomPairing() {
        return atomPairing;
    }

    @Override
    public Transformation getTransformation() {
        return transformation;
    }

    @Override
    public double getAlignmentScore() {
        return alignmentScore;
    }
}
