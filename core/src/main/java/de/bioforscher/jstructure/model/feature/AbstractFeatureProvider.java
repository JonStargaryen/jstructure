package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;

/**
 * The abstract implementation of each {@link FeatureProvider}. Implicitly registers each provider in a global
 * "registry" implemented by {@link FeatureProviderRegistry}.
 * Created by S on 16.01.2017.
 */
public abstract class AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(AbstractFeatureProvider.class);

    /**
     * Runs the computations implemented by this feature provider and assigns the results to the given container.
     * @param protein the container to process
     */
    public void process(Protein protein) {
        FeatureProvider annotation = getClass().getDeclaredAnnotation(FeatureProvider.class);
        Group firstGroup = protein.getGroups().get(0);
        // additional features must be computed beforehand
        if (!Arrays.equals(annotation.requiredFeatures(), new String[0])) {
            for(String requiredFeature : annotation.requiredFeatures()) {
                // feature present
                //TODO maybe move back to some feature container root
                if(firstGroup.getFeatureMap().containsKey(requiredFeature)) {
                    continue;
                }

                // need to compute
                AbstractFeatureProvider resolvedFeatureProvider = FeatureProviderRegistry.getInstance().resolve(requiredFeature);
                logger.info("computing {}: using {} to compute required feature {}",
                        Arrays.toString(annotation.providedFeatures()),
                        resolvedFeatureProvider.getClass().getSimpleName(),
                        requiredFeature);
                resolvedFeatureProvider.process(protein);
            }
        }

        processInternally(protein);

        //TODO define some clean-up routine?
    }

    protected abstract void processInternally(Protein protein);
}
