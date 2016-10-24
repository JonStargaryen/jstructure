package de.bioforscher.jstructure.model.structure;

import java.util.stream.Stream;

/**
 * Specifies the capabilities of a chain container (mostly a {@link Protein}).
 * Created by S on 30.09.2016.
 */
public interface ChainContainer extends GroupContainer {
    /**
     * Access to all chain objects associated to this container.
     * @return a stream of chains
     */
    Stream<Chain> chains();

    //TODO functionality to access aminoAcid chains, ligands, nucleotides etc
    //TODO chainFilter class
}
