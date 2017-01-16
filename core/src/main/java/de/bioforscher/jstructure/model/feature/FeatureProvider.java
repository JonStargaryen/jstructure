package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.Protein;

import java.util.List;

/**
 * Specification of an algorithm which can compute features and write them to the feature map.
 * Created by S on 02.10.2016.
 */
public interface FeatureProvider {
    /**
     *
     * @param protein
     */
    void process(Protein protein);

    /**
     *
     * @return
     */
    List<String> getProvidedFeaturesNames();


    //TODO define some clean-up routine
}