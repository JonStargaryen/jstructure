package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.Identifiable;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;

/**
 * Created by bittrich on 6/19/17.
 */
public class AbstractBasicContainer extends AbstractFeatureable implements OriginIndicator, Identifiable {
    private final Protein origin;
    private String identifier;

    AbstractBasicContainer(Protein origin) {
        this.origin = origin;
        this.identifier = origin.getIdentifier();
    }

    AbstractBasicContainer(AbstractBasicContainer container) {
        this.origin = container.origin;
        this.identifier = container.identifier;
    }

    public Protein getOrigin() {
        return origin;
    }

    @Override
    public String getIdentifier() {
        return identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}