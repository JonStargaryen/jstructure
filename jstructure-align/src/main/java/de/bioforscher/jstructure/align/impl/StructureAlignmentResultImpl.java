package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.StructureAlignmentResult;
import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * The result of an structure alignment.
 * Created by bittrich on 6/19/17.
 */
public class StructureAlignmentResultImpl implements StructureAlignmentResult {
    private final GroupContainer originalReference;
    private final GroupContainer originalQuery;
    private final GroupContainer alignedQuery;
    private final Transformation transformation;
    private final double alignmentScore;

    StructureAlignmentResultImpl(GroupContainer originalReference,
                                 GroupContainer originalQuery,
                                 GroupContainer alignedQuery,
                                 Transformation transformation,
                                 double alignmentScore) {
        this.originalReference = originalReference;
        this.originalQuery = originalQuery;
        this.alignedQuery = alignedQuery;
        this.transformation = transformation;
        this.alignmentScore = alignmentScore;
    }

    @Override
    public GroupContainer getOriginalReference() {
        return originalReference;
    }

    @Override
    public GroupContainer getOriginalQuery() {
        return originalQuery;
    }

    @Override
    public GroupContainer getAlignedQuery() {
        return alignedQuery;
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
