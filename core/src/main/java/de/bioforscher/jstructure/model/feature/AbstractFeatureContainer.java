package de.bioforscher.jstructure.model.feature;

import de.bioforscher.jstructure.model.structure.container.StructureContainer;

import java.util.HashMap;
import java.util.Map;

/**
 * The abstract implementation of {@link StructureContainer}. Especially providing access to the featureMap and standardizing its
 * usage.
 * Created by S on 15.11.2016.
 */
public abstract class AbstractFeatureContainer implements StructureContainer {
    private Map<Enum, Object> featureMap;

    public AbstractFeatureContainer() {
        this.featureMap = new HashMap<>();
    }

    //TODO maybe some unique name for each impl - protein:pdbName, chain:chainId, group:pdbName-resNum, atom:element-pdbSerial

    @Override
    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}
