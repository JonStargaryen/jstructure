package de.bioforscher.jstructure.align;

import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.List;

/**
 * The result of a structural alignment.
 * Created by S on 12.07.2017.
 */
public interface StructureAlignmentResult {
    /**
     * The group container used as reference in this alignment.
     * @return the original reference as {@link AtomContainer}
     */
    AtomContainer getOriginalReference();

    /**
     * The group container which was supposed to be aligned to the reference.
     * @return the original query as {@link AtomContainer}
     */
    AtomContainer getOriginalQuery();

    /**
     *
     * @return
     */
    AtomContainer getAlignedReference();

    /**
     * The aligned query container - may be only a subset of entities from the original query.
     * @return a manipulated {@link AtomContainer}
     */
    AtomContainer getAlignedQuery();

    List<Pair<Atom, Atom>> getAtomPairing();

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
