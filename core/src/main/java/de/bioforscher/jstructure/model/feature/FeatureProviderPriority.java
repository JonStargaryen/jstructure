package de.bioforscher.jstructure.model.feature;

/**
 * The annotation used to create {@link FeatureProvider} and track them by the {@link FeatureServiceRegistry}.
 * Created by S on 16.01.2017.
 */
@interface FeatureProviderPriority {
    /**
     * The default priority of each and every {@link FeatureProvider}.
     */
    int DEFAULT_PRIORITY = 100;

    /**
     * The priority of how this {@link FeatureProvider} is utilized to compute requested features. Lower values
     * indicate services of higher priority which will be used instead of providers of lower priority. When the priority
     * is tied, the resolved {@link FeatureProvider} is random.
     * @return the priority of this service
     */
    int priority() default DEFAULT_PRIORITY;
}
