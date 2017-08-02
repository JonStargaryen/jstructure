package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Structure;

import java.io.*;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The abstract implementation of each {@link FeatureProvider}. Feature provider implementations are expected to be
 * stateless.
 * TODO a register function which fills each featureMap with the provided keys of this provider, thus multiple predictors could be employed concurrently
 * TODO strict mode flag
 * TODO indicate where the results are attached to - on atom level? residues? chains? whole protein? a mix?
 *  Created by S on 16.01.2017.
 */
public abstract class FeatureProvider {
    /**
     * Runs the computations implemented by this feature provider and assigns the results to the given container.
     * @param structure the container to process
     */
    public void process(Structure structure) {
        preprocessInternally(structure);

        // delegate to concrete implementation
        processInternally(structure);

        // 'hook' to postprocessing routine, empty by default but can be overridden by implementations if needed
        postprocessInternally(structure);

        // registered the computed feature entry class
        provides().forEach(structure::registerFeature);
    }

    protected List<Class<? extends FeatureContainerEntry>> provides() {
        return Collections.emptyList();
    }

    private void preprocessInternally(Structure structure) {

    }

    /**
     * Some postprocess method with will be executed on the container after the computation has finished. Override when
     * there is a need to clean-up some variables etc in the container.
     * @param structure the container to clean
     */
    protected void postprocessInternally(Structure structure) {

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

    protected abstract void processInternally(Structure structure);
}
