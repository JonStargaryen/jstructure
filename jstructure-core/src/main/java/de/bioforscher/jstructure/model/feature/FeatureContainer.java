package de.bioforscher.jstructure.model.feature;

import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The internal feature map class associated to instances implementing {@link Featureable}.
 * Created by S on 28.04.2017.
 */
public class FeatureContainer {
    private List<FeatureContainerEntry> features;

    public FeatureContainer() {
        this.features = new ArrayList<>();
    }

    public void addFeature(FeatureContainerEntry entry) {
        features.add(entry);
    }

    /**
     * Access particular features, identified by its content type.
     * @param contentClass the class of the content of interest
     * @param <C>
     * @throws NoSuchElementException occurs when no entry is present
     * @return the <b>first</b> relevant
     */
    public <C extends FeatureContainerEntry> C getFeature(Class<C> contentClass) {
        return getFeatureOptional(contentClass)
                .orElseThrow(() -> new NoSuchElementException("no feature entry for content class '" +
                        contentClass.getSimpleName()));
    }

    /**
     * The safe way to access particular features, identified by its entry type.
     * @param contentClass the class of the content of interest
     * @param <C>
     * @return the <b>first</b> relevant, wrapped as {@link Optional}
     */
    public <C extends FeatureContainerEntry> Optional<C> getFeatureOptional(Class<C> contentClass) {
        return filterByContent(contentClass)
                .findFirst();
    }

    /**
     * Returns all associated entries of a given type.
     * @param contentClass the class of the content of interest
     * @return a collection of relevant entries
     */
    public <C extends FeatureContainerEntry> List<C> getFeatures(Class<C> contentClass) {
        return filterByContent(contentClass)
                .collect(Collectors.toList());
    }

    private <C extends FeatureContainerEntry> Stream<C> filterByContent(Class<C> contentClass) {
        return features.stream()
                .filter(contentClass::isInstance)
                .map(contentClass::cast);
    }

    /**
     * Low-level access to the entries.
     * @return the delegated list
     */
    public List<FeatureContainerEntry> getFeatures() {
        return features;
    }

    /**
     * Set the internal list to a specific value. Used after cloning of model instances.
     * @param features the new value of the feature container
     */
    public void setFeatures(List<FeatureContainerEntry> features) {
        this.features = features;
    }
}
