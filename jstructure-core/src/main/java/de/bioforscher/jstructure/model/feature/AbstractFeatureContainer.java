package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.container.StructureContainer;

import java.util.HashMap;
import java.util.Map;

/**
 * The abstract implementation of {@link StructureContainer}. Especially providing access to the featureMap and
 * standardizing its usage.
 * Created by S on 15.11.2016.
 */
public abstract class AbstractFeatureContainer implements StructureContainer {
    private Map<String, Object> featureMap;

    protected AbstractFeatureContainer() {
        this.featureMap = new HashMap<>();
    }

    @Override
    public Map<String, Object> getFeatureMap() {
        return featureMap;
    }

    @Override
    public void setFeatureMap(Map<String, Object> map) {
        this.featureMap = map;
    }
}
