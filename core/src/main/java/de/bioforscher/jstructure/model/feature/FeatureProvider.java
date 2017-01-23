package de.bioforscher.jstructure.model.feature;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Annotation of an algorithm which can compute features and write them to the feature map of a {@link FeatureContainer}.
 * TODO distinction between global and local features?
 * Created by S on 02.10.2016.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface FeatureProvider {
    int DEFAULT_PRIORITY = 100;

    FeatureType type() default FeatureType.ANNOTATION;

    /**
     * The priority of this feature provider which can be used to assign an ordering on which implementation is employed
     * to compute certain features. Override this method to assign custom priorities.
     * @return this feature providers priority
     */
    int priority() default DEFAULT_PRIORITY;

    /**
     * All features provided by the implementing class. This list should contain at least 1 element.
     * @return all features provided by this implemented
     */
    String[] provides();

    /**
     * Access to all required features which need to be present in order for the computation about to happen to succeed.
     * E.g., the computation of something like the loop fraction requires the previous calculation/annotation of all
     * secondary structure information. By resolving these dependencies, the {@link AbstractFeatureProvider} will try to
     * compute all requirements beforehand automatically. This list can be empty.
     * @return all features needed for the computation implemented by this provider
     */
    String[] requires() default {};
}