package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.Protein;

/**
 * Created by bittrich on 6/19/17.
 */
interface OriginIndicator {
    /**
     * The model instance on which this container was created upon.
     * @return access to the origin of this container
     */
    Protein getOrigin();
}