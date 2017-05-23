package de.bioforscher.jstructure.model.structure.selection;

import java.util.NoSuchElementException;

/**
 * Thrown when a selection cannot find any suitable element but was forced to do so by a call to
 * {@link Selection.AtomSelection#asAtom()}, {@link Selection.GroupSelection#asGroup()} or
 * {@link Selection.ChainSelection#asChain()}. *
 * Created by bittrich on 5/23/17.
 */
public class SelectionException extends NoSuchElementException {
    public SelectionException(String s) {
        super(s);
    }
}
