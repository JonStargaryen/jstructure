package de.bioforscher.jstructure.model.feature;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Annotation of an algorithm which can calculate features and write them to the feature map of a {@link Featureable}.
 * Created by S on 02.10.2016.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface FeatureProvider {
    int DEFAULT_PRIORITY = 100;

    FeatureOrigin origin() default FeatureOrigin.ANNOTATION;

    /**
     * The priority of this feature provider which can be used to assign an ordering on which implementation is employed
     * to calculate certain features. Override this method to assign custom priorities.
     * @return this feature providers priority
     */
    int priority() default DEFAULT_PRIORITY;

    /**
     * All features provided by the implementing class. This list should contain at least 1 element.
     * @return all features provided by this implemented
     */
    Class<? extends FeatureContainerEntry>[] provides();

    /**
     * Access to all required features which need to be present in order for the computation about to happen to succeed.
     * E.g., the computation of something like the loop fraction requires the previous calculation/annotation of all
     * secondary structure information. By resolving these dependencies, the {@link AbstractFeatureProvider} will try to
     * calculate all requirements beforehand automatically. This list can be empty.
     * @return all features needed for the computation implemented by this provider
     */
    Class<? extends FeatureContainerEntry>[] requires() default {};

    /**
     * The possible feature types - predictions and annotations.
     */
    enum FeatureOrigin {
        ANNOTATION,
        PREDICTION
    }
}