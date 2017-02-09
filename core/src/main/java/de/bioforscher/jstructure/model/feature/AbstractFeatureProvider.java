package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Objects;
import java.util.stream.Stream;

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
        if (!Arrays.equals(annotation.requires(), new String[0])) {
            for(String requiredFeature : annotation.requires()) {
                // feature present
                //TODO maybe move back to some feature container root
                if(firstGroup.getFeatureMap().containsKey(requiredFeature)) {
                    continue;
                }

                // need to compute
                AbstractFeatureProvider resolvedFeatureProvider = FeatureProviderRegistry.resolve(requiredFeature);
                logger.info("computing {}: using {} to compute required feature {}",
                        Arrays.toString(annotation.provides()),
                        resolvedFeatureProvider.getClass().getSimpleName(),
                        requiredFeature);
                resolvedFeatureProvider.process(protein);
            }
        }

        // delegate to concrete implementation
        processInternally(protein);

        // 'hook' to postprocessing routine, empty by default but can be overridden by implementations if needed
        postprocessInternally(protein);
    }

    /**
     * Some postprocess method with will be executed on the container after the computation has finished. Override when
     * there is a need to clean-up some variables etc in the container.
     * @param protein the container to clean
     */
    protected void postprocessInternally(Protein protein) {

    }

    protected InputStream getResourceAsStream(String filepath) {
        return Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResourceAsStream(filepath),
                "failed to find resource as InputStream");
    }

    protected Stream<String> getLinesFromResource(String filepath) {
        return new BufferedReader(new InputStreamReader(getResourceAsStream(filepath))).lines();
    }

    protected abstract void processInternally(Protein protein);
}
