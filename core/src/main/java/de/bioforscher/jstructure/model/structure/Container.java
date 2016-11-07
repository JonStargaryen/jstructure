package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.feature.FeatureContainer;
import de.bioforscher.jstructure.feature.FeatureProvider;

/**
 * Marks this object as part of a hierarchy. Specifies no methods directly, but is useful so
 * {@link FeatureProvider} can be generically tied to classes implementing this
 * interface.
 * Created by S on 02.10.2016.
 */
public interface Container extends FeatureContainer {
    //TODO maybe we should drop clone all together and provide a constructor such as Cloneable(Cloneable instance)
//    /**
//     * Clones this container object. <b>Caution:</b> the <code>featureMap</code> as well as the parent container is
//     * shared with the original instance.
//     * @return a deep copy of this object, shallow copy of the parent container (if any) and the <code>featureMap</code>
//     */
//    Container clone();
}
