package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

/**
 * Thrown when a selection cannot find any suitable element but was forced to do so by a call to
 * {@link Selection.AtomSelection#asAtom()}, {@link Selection.GroupSelection#asGroup()} or
 * {@link Selection.ChainSelection#asChain()}. Also happening when {@link Selection} was asked to perform operations on
 * selected elements which were actually not compatible with the selection, i.e. casting normal groups to
 * {@link AminoAcid}.
 * Created by bittrich on 5/23/17.
 */
public class SelectionException extends RuntimeException {
    public SelectionException(String s) {
        super(s);
    }

    public SelectionException(String message, Throwable cause) {
        super(message, cause);
    }

    public SelectionException(Throwable cause) {
        super(cause);
    }

    public SelectionException() {
    }
}
