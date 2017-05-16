package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.structure.container.StructureContainer;

/**
 * Specifies the capability of a container to return a specific name. This is not necessarily corresponding to data
 * parsed from <code>PDB</code> files.
 * Created by S on 18.11.2016.
 */
public interface Identifiable {
    /**
     * A more-or-less unique identifier. Specifically, this make {@link Object#toString()} calls more informative on
     * {@link StructureContainer} objects created by selections which do not
     * necessarily feature a identifier parsed from a <code>PDB</code> file.
     * @return the source of this container
     */
    String getIdentifier();

    void setIdentifier(String identifier);
}
