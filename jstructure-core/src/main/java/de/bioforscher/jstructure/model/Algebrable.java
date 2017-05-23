package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;

/**
 * The implementing class supports algebraic operations.
 * Created by bittrich on 5/23/17.
 */
public interface Algebrable<T extends LinearAlgebra.AlgebraicOperations> {
    /**
     * Access to algebraic operations on this instance.
     * @return a {@link de.bioforscher.jstructure.mathematics.LinearAlgebra} instance
     */
    T algebra();
}
