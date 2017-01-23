package feature;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfilePredictor;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import org.junit.Assert;
import org.junit.Test;

import java.util.NoSuchElementException;

/**
 * The functional test for the feature service registry.
 * Created by bittrich on 1/17/17.
 */
public class FeatureProviderRegistryTest {
    private static final String NOT_REGISTERED_FEATURE_NAME = "NOT_REGISTERED_FEATURE_NAME";

    @Test
    public void shouldGetAllRegisteredFeatureProviders() throws Exception {
        FeatureProviderRegistry.getRegisteredFeatureProviders().forEach(System.out::println);
    }

    @Test
    public void shouldGetSupportedFeatures() throws Exception {
        FeatureProviderRegistry.getSupportedFeatures().forEach(System.out::println);
    }

    @Test(expected = NoSuchElementException.class)
    public void shouldFailToResolve() throws Exception {
        FeatureProviderRegistry.resolve(NOT_REGISTERED_FEATURE_NAME);
    }

    @Test
    public void shouldResolveProvider() {
        //TODO fails in the core-module because other modules are not visible - can this be actually fixed?
        FeatureProviderRegistry.resolve("MEMBRANE");
    }

    @Test
    public void shouldResolvePredictor() {
        Assert.assertTrue(FeatureProviderRegistry.resolvePredictor(EnergyProfileCalculator.SOLVATION_ENERGY).getClass().equals(EnergyProfilePredictor.class));
    }
}