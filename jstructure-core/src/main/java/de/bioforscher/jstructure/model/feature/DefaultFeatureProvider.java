package de.bioforscher.jstructure.model.feature;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface DefaultFeatureProvider {
    Class<? extends FeatureProvider> value();
}
