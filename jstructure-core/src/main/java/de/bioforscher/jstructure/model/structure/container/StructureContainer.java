package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.Copyable;
import de.bioforscher.jstructure.model.feature.*;
import de.bioforscher.jstructure.model.structure.AtomRecordWriter;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;


/**
 * Marks this object as part of a hierarchy. Specifies no methods directly, but is useful so
 * {@link FeatureProvider} can be generically tied to classes implementing this
 * interface. As a rule of thumb, there are concrete groupings of elements (such as {@link Group} or {@link Chain} and
 * rather loose collections of similar elements such as {@link AtomContainer} which does not necessarily assumes any
 * clearAtoms connection between a set of atoms.
 * Created by S on 02.10.2016.
 */
public interface StructureContainer extends Featureable, AtomRecordWriter, Copyable {
    Logger logger = LoggerFactory.getLogger(StructureContainer.class);

    String getIdentifier();

    void setIdentifier(String identifier);

    Structure getParentStructure();

    /**
     * Access particular features, identified by its content type.
     * @param contentClass the class of the content of interest
     * @param <C> the content class returned
     * @return the <b>first</b> relevant
     */
    @Override
    default <C extends FeatureContainerEntry> C getFeature(Class<C> contentClass) {
        Optional<C> optional = getFeatureContainer().getFeatureOptional(contentClass);
        if(optional.isPresent()) {
            return optional.get();
        } else {
            FeatureProvider featureProvider = DefaultFeatureProviderMap.resolve(contentClass);
            logger.debug("feature {} was not present, using {} to compute",
                    contentClass.getSimpleName(),
                    featureProvider.getClass().getSimpleName());
            featureProvider.process(getParentStructure());
            return getFeatureContainer().getFeatureOptional(contentClass).orElseThrow(() ->
                    new ComputationException("feature " + contentClass.getSimpleName() + " could not be computed by " +
                            featureProvider.getClass().getSimpleName()));
        }
    }

    /**
     * Handles the resolution of feature providers. Contains a map of initialized feature provider instances which can,
     * thus, be reused and must be stateless.
     */
    class DefaultFeatureProviderMap {
        private static final Logger logger = LoggerFactory.getLogger(DefaultFeatureProviderMap.class);
        private static final Map<Class<? extends FeatureContainerEntry>, FeatureProvider> featureProviderMap = new HashMap<>();

        static FeatureProvider resolve(Class<? extends FeatureContainerEntry> featureContainerEntry) {
            if(!featureProviderMap.containsKey(featureContainerEntry)) {
                try {
                    Class<? extends FeatureProvider> featureProviderClass = featureContainerEntry.getAnnotation(DefaultFeatureProvider.class).value();
                    FeatureProvider featureProviderInstance = featureProviderClass.newInstance();
                    featureProviderMap.put(featureContainerEntry, featureProviderInstance);
                    logger.debug("establishing mapping {} => {}",
                            featureContainerEntry.getSimpleName(),
                            featureProviderInstance.getClass().getSimpleName());
                    return featureProviderInstance;
                } catch (NullPointerException e) {
                    throw new ComputationException("missing DefaultFeatureProvider annotation for class " +
                            featureContainerEntry.getSimpleName() + " - cannot resolve feature provider");
                } catch (Exception e) {
                    throw new ComputationException(e);
                }
            }
            return featureProviderMap.get(featureContainerEntry);
        }
    }
}
