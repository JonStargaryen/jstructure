package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.List;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a getChain container (mostly a {@link Protein}).
 * Created by S on 30.09.2016.
 */
public interface ChainContainer extends GroupContainer {
    /**
     * Never manipulate the returned collection as it is not guaranteed the actually modify the internal list(s).
     * @return all associated chains
     */
    List<Chain> getChains();

    /**
     * Access to all getChain objects associated to this container.
     * @return a select of chains
     */
    default Stream<Chain> chains() {
        return getChains().stream();
    }

    default ChainContainer getCopy() {
        return (ChainContainer) GroupContainer.super.getCopy();
    }
}
