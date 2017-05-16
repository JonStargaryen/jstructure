package de.bioforscher.jstructure.model.feature;

import org.reflections.Reflections;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The global registry of known {@link FeatureProvider} instances. Every {@link FeatureProvider} automatically registers
 * itself in this list, so the user can retrieve the default implementation for certain features (identified by name)
 * and also for a {@link FeatureProvider} to resolve other implementations when it depends on them without explicitly
 * referencing them or directly depending on them or the concrete implementation (as many features can be computed and
 * defined in a various number of ways).
 * Created by S on 16.01.2017.
 */
public class FeatureProviderRegistry {
    private static final Logger logger = LoggerFactory.getLogger(FeatureProviderRegistry.class);
    private static final FeatureProviderRegistry INSTANCE = new FeatureProviderRegistry();
    private Map<String, TreeMap<Integer, AbstractFeatureProvider>> registeredFeatureProviders;

    private FeatureProviderRegistry() {
        logger.info("setting up feature provider registry");
        registeredFeatureProviders = new HashMap<>();

        scanForFeatureProviders();
    }

    private void scanForFeatureProviders() {
        Reflections reflections = new Reflections("de.bioforscher");
        Set<Class<?>> annotatedFeatureProviders = reflections.getTypesAnnotatedWith(FeatureProvider.class);

        for(Class<?> annotatedFeatureProvider : annotatedFeatureProviders) {
            FeatureProvider annotation = annotatedFeatureProvider.getDeclaredAnnotation(FeatureProvider.class);
            logger.info("registering provider {} with priority {}",
                    annotatedFeatureProvider.getSimpleName(),
                    annotation.priority());

            for(String providedFeature : annotation.provides()) {
                TreeMap<Integer, AbstractFeatureProvider> providers = registeredFeatureProviders.getOrDefault(providedFeature, new TreeMap<>());
                try {
                    int priority = annotation.priority();
                    if(providers.containsKey(priority)) {
                        providers.put(priority, (AbstractFeatureProvider) annotatedFeatureProvider.newInstance());
                    } else {
                        // decrease priority to dodge duplicates, as only 1 would be registered
                        while(providers.containsKey(priority)) {
                            priority++;
                        }

                        providers.put(priority, (AbstractFeatureProvider) annotatedFeatureProvider.newInstance());
                    }
                } catch (InstantiationException | IllegalAccessException e) {
                    throw new InstantiationError("could not create instance of feature provider " + annotatedFeatureProvider.getSimpleName());
                }
                registeredFeatureProviders.put(providedFeature, providers);
            }
        }
    }

    /**
     * A lists all registered {@link FeatureProvider}.
     * @return all registered services
     */
    public static List<AbstractFeatureProvider> getRegisteredFeatureProviders() {
        return INSTANCE.registeredFeatureProviders.entrySet().stream()
                .map(Map.Entry::getValue)
                .map(Map::entrySet)
                .flatMap(Collection::stream)
                .map(Map.Entry::getValue)
                .distinct()
                .collect(Collectors.toList());
    }

    /**
     * Access to all features supported by the current setup.
     * @return a list of all features names for which at least 1 {@link FeatureProvider} is registered which can compute
     * it
     */
    public static List<String> getSupportedFeatures() {
        return INSTANCE.registeredFeatureProviders.entrySet().stream()
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
    }

    //TODO builder-esque resolving

    public static AbstractFeatureProvider resolveAnnotator(String requestedFeature) {
        return resolveInternal(requestedFeature)
                .filter(provider -> provider.getClass().getAnnotation(FeatureProvider.class).type().equals(FeatureType.ANNOTATION))
                .findFirst()
                .orElseThrow(() -> noProvider(requestedFeature));
    }

    public static AbstractFeatureProvider resolvePredictor(String requestedFeature) {
        return resolveInternal(requestedFeature)
                .filter(provider -> provider.getClass().getAnnotation(FeatureProvider.class).type().equals(FeatureType.PREDICTION))
                .findFirst()
                .orElseThrow(() -> noProvider(requestedFeature));
    }

    private static Stream<AbstractFeatureProvider> resolveInternal(String requestedFeature) {
        try {
            return INSTANCE.registeredFeatureProviders.get(requestedFeature).entrySet().stream().map(Map.Entry::getValue);
        } catch (NullPointerException e) {
            throw noProvider(requestedFeature);
        }
    }

    /**
     * Allows to resolve arbitrary feature name and returns the {@link FeatureProvider} which should be employed to
     * compute it.
     * @return the favored {@link FeatureProvider} (according to the priority with which that feature provider is
     * registered)
     * @throws java.util.NoSuchElementException if no {@link FeatureProvider} is registered which can compute the
     * requested feature
     */
    public static AbstractFeatureProvider resolve(String requestedFeature) {
        return resolveInternal(requestedFeature).findFirst().orElseThrow(() -> noProvider(requestedFeature));
    }

    private static NoSuchElementException noProvider(String requestedFeature) {
        throw new NoSuchElementException("no provider is registered for '" + requestedFeature + "'");
    }
}