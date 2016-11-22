package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Protein;

/**
 * Specification of an algorithm which can compute features and write them to the feature map.
 * Created by S on 02.10.2016.
 */
public interface FeatureProvider {
    void process(Protein protein);

    //TODO define some clean-up routine
    //TODO standardize the featureNames they provide
    //TODO what features does this algorithm depend on?
}