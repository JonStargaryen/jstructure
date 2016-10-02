package de.bioforscher.jstructure.feature;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;

/**
 *
 * @param <F> the minimal container this feature asa is applicable for - all parent containers up to the root can
 *           be handled, but none of the potential child container (e.g. a FeatureProvider<{@link Chain}> is able of
 *           processing a {@link Protein}, but no {@link Residue}
 * Created by S on 02.10.2016.
 */
public interface FeatureProvider<F extends FeatureContainer> {
    void process(F featureContainer);
}