package de.bioforscher.jstructure.model.structure;

import java.util.NoSuchElementException;
import java.util.Optional;
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

    /**
     *
     * @param chainId
     * @return
     */
    default Optional<Chain> chain(final String chainId) {
        return chains().filter(chain -> chain.getChainId().equals(chainId)).findFirst();
    }

    /**
     *
     * @param chainId
     * @param residueNumber
     * @return
     */
    default Optional<Residue> residue(final String chainId, final int residueNumber) {
        return chain(chainId).orElseThrow(NoSuchElementException::new).residue(residueNumber);
    }
}
