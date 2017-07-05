package de.bioforscher.jstructure.feature;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfilePredictor;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import org.junit.Assert;
import org.junit.Test;

import java.util.NoSuchElementException;

/**
 * The functional test for the feature service registry.
 * Created by bittrich on 1/17/17.
 */
public class FeatureProviderRegistryTest {
    class NotRegisteredFeature extends FeatureContainerEntry {
        public NotRegisteredFeature(AbstractFeatureProvider featureProvider) {
            super(featureProvider);
        }
    }

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
        FeatureProviderRegistry.resolve(NotRegisteredFeature.class);
    }

    @Test
    public void shouldResolveProvider() {
        FeatureProviderRegistry.resolve(AccessibleSurfaceArea.class);
    }

    @Test
    public void shouldResolvePredictor() {
        Assert.assertTrue(FeatureProviderRegistry.resolvePredictor(EnergyProfile.class).getClass().equals(EnergyProfilePredictor.class));
    }
}