package feature;

import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import org.junit.Test;

import java.util.NoSuchElementException;

/**
 * The functional test for the feature service registry.
 * Created by bittrich on 1/17/17.
 */
public class FeatureProviderRegistryTest {
    private FeatureProviderRegistry instance = FeatureProviderRegistry.getInstance();
    private static final String NOT_REGISTERED_FEATURE_NAME = "NOT_REGISTERED_FEATURE_NAME";

    @Test
    public void shouldGetAllRegisteredFeatureProviders() throws Exception {
        instance.getRegisteredFeatureProviders().forEach(System.out::println);
    }

    @Test
    public void shouldGetSupportedFeatures() throws Exception {
        instance.getSupportedFeatures().forEach(System.out::println);
    }

    @Test(expected = NoSuchElementException.class)
    public void shouldFailToResolve() throws Exception {
        instance.resolve(NOT_REGISTERED_FEATURE_NAME);
    }

    @Test
    public void shouldResolveProvider() {
        instance.resolve();
    }
}