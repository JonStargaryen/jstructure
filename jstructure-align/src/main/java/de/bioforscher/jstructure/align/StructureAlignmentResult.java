package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

/**
 * The result of a structural alignment.
 * Created by S on 12.07.2017.
 */
public interface StructureAlignmentResult {
    /**
     * The group container used as reference in this alignment.
     * @return the original reference as {@link GroupContainer}
     */
    GroupContainer getOriginalReference();

    /**
     * The group container which was supposed to be aligned to the reference.
     * @return the original query as {@link GroupContainer}
     */
    GroupContainer getOriginalQuery();

    /**
     * The aligned query container - may be only a subset of entities from the original query.
     * @return a manipulated {@link GroupContainer}
     */
    GroupContainer getAlignedQuery();

    /**
     * The transformation which can recreate this alignment.
     * @return a {@link Transformation} instance
     */
    Transformation getTransformation();

    /**
     * The root-mean square deviation which this alignment achieved.
     * @return a double score
     */
    double getAlignmentScore();
}
