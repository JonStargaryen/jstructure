package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;

import java.util.List;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a getChain container (mostly a {@link Structure}).
 * Created by S on 30.09.2016.
 */
public interface ChainContainer extends GroupContainer {
    /**
     * Never manipulate the returned collection as it is not guaranteed the actually modify the internal list(s).
     * @return all associated chains
     */
    List<Chain> getChains();

    /**
     * Access to all chain objects associated to this container.
     * @return a stream of chains
     */
    default Stream<Chain> chains() {
        return getChains().stream();
    }

    /**
     * Access to all chains which actually contains amino acids.
     * @return a stream of chains with at least 1 amino acid
     */
    default Stream<Chain> chainsWithAminoAcids() {
        return chains()
                .filter(chain -> chain.aminoAcids().count() > 0);
    }
}
