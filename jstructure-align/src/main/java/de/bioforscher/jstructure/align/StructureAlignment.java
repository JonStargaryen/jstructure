package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * The result of an alignment.
 * Created by bittrich on 6/19/17.
 */
public class StructureAlignment {
    private final GroupContainer originalReference;
    private final GroupContainer originalQuery;
    private final GroupContainer alignedQuery;
    private final Transformation transformation;
    private final double alignmentScore;

    StructureAlignment(GroupContainer originalReference,
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

    public GroupContainer getOriginalReference() {
        return originalReference;
    }

    public GroupContainer getOriginalQuery() {
        return originalQuery;
    }

    public GroupContainer getAlignedQuery() {
        return alignedQuery;
    }

    public Transformation getTransformation() {
        return transformation;
    }

    public double getAlignmentScore() {
        return alignmentScore;
    }
}
