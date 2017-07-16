package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The abstract implementation of each {@link FeatureProvider}. Implicitly registers each provider in a global
 * "registry" implemented by {@link FeatureProviderRegistry}. Feature provider implementations are expected to be
 * stateless.
 * TODO a register function which fills each featureMap with the provided keys of this provider, thus multiple predictors could be employed concurrently
 * TODO strict mode flag
 * TODO indicate where the results are attached to - on atom level? residues? chains? whole protein? a mix?
 *  Created by S on 16.01.2017.
 */
public abstract class AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(AbstractFeatureProvider.class);

    protected AbstractFeatureProvider() {
        // register feature providers upon creation in the registry - TODO does not fix as classes are not guaranteed to be loaded at all
//        FeatureProviderRegistry.register(this.getClass());
    }

    /**
     * Runs the computations implemented by this feature provider and assigns the results to the given container.
     * @param protein the container to processUniProtId
     */
    public void process(Structure protein) {
        FeatureProvider annotation = getClass().getDeclaredAnnotation(FeatureProvider.class);
        // additional features must be computed beforehand
        if (!Arrays.equals(annotation.requires(), new String[0])) {
            for(Class<? extends FeatureContainerEntry> requiredFeature : annotation.requires()) {
                // feature present
                if(protein.featureIsPresent(requiredFeature)) {
                    continue;
                }

                // need to calculate
                AbstractFeatureProvider resolvedFeatureProvider = FeatureProviderRegistry.resolve(requiredFeature);
                logger.debug("computing {}: using {} to calculate required feature {}",
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

        // registered the computed feature entry class
        Stream.of(annotation.provides())
                .forEach(protein::registerFeature);
    }

    /**
     * Some postprocess method with will be executed on the container after the computation has finished. Override when
     * there is a need to clean-up some variables etc in the container.
     * @param protein the container to clean
     */
    protected void postprocessInternally(Structure protein) {

    }

    protected static InputStream getResourceAsInputStream(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        InputStream is = ccl.getResourceAsStream(filename);
        return Objects.requireNonNull(is);
    }

    protected static Stream<String> getResourceAsStream(String filename) {
        return getResourceAsLines(filename).stream();
    }

    protected static List<String> getResourceAsLines(String filename) {
        try {
            try (InputStreamReader inputStreamReader = new InputStreamReader(getResourceAsInputStream(filename))) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    return bufferedReader.lines()
                            .collect(Collectors.toList());
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    protected abstract void processInternally(Structure protein);
}
