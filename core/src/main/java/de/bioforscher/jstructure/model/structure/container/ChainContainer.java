package de.bioforscher.jstructure.model.structure.container;

import com.sun.org.apache.regexp.internal.RE;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a getChain container (mostly a {@link Protein}).
 * Created by S on 30.09.2016.
 */
public interface ChainContainer extends GroupContainer {
    List<Chain> getChains();

    /**
     * Access to all getChain objects associated to this container.
     * @return a stream of chains
     */
    default Stream<Chain> chains() {
        return getChains().stream();
    }

    /**
     *
     * @param chainId
     * @return
     */
    default Chain getChain(final String chainId) {
        return findChain(chainId).orElseThrow(NoSuchElementException::new);
    }

    /**
     *
     * @param chainId
     * @param residueNumber
     * @return
     */
    default Residue getResidue(final String chainId, final int residueNumber) {
        return findResidue(chainId, residueNumber).orElseThrow(NoSuchElementException::new);
    }

    default Optional<Residue> findResidue(final String chainId, final int residueNumber) {
        return findChain(chainId).orElseThrow(NoSuchElementException::new).findResidue(residueNumber);
    }

    default Optional<Chain> findChain(String chainId) {
        return chains().filter(chain -> chain.getChainId().equals(chainId)).findFirst();
    }
}
