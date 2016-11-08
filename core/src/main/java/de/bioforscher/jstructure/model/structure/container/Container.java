package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.feature.FeatureContainer;
import de.bioforscher.jstructure.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.AtomRecordWriter;

/**
 * Marks this object as part of a hierarchy. Specifies no methods directly, but is useful so
 * {@link FeatureProvider} can be generically tied to classes implementing this
 * interface.
 * Created by S on 02.10.2016.
 */
public interface Container extends FeatureContainer, AtomRecordWriter {

}
