package de.bioforscher.jstructure.model.feature;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The global registry of known {@link FeatureProvider} instances. Every {@link FeatureProvider} automatically registers
 * itself in this list, so the user can retrieve the default implementation for certain features (identified by name)
 * and also for a {@link FeatureProvider} to resolve other implementations when it depends on them without explicitly
 * referencing them or directly depending on them or the concrete implementation (as many features can be computed and
 * defined in a various number of ways).
 * Created by S on 16.01.2017.
 */
public class FeatureServiceRegistry {
    private static final FeatureServiceRegistry INSTANCE = new FeatureServiceRegistry();
    /**
     *
     */
    private Map<String, Map<Integer, FeatureProvider>> registeredFeatureProviders;

    private FeatureServiceRegistry() {
        this.registeredFeatureProviders = new ConcurrentHashMap<>();
    }

    /**
     * Access to this classes
     * @return the singleton instance of this class
     */
    public static FeatureServiceRegistry getInstance() {
        return INSTANCE;
    }

    public Map<String, Map<Integer, FeatureProvider>> getFeatureProviderMap() {
        return registeredFeatureProviders;
    }

    /**
     * A lists all registered {@link FeatureProvider}.
     * @return all registered services
     */
    public List<FeatureProvider> getAllRegisteredFeatureProviders() {
        return null;
    }

    /**
     * Registers a given {@link FeatureProvider}.
     * @param featureProviderToRegister the instance to add
     */
    public void registerFeatureProvider(FeatureProvider featureProviderToRegister) {
//        registeredFeatureProviders.add(featureProviderToRegister);
    }
}
