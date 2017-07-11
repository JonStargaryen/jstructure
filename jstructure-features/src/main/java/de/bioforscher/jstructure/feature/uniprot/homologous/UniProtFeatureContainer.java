package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Information obtained by a UniProt mapping to homologous protein chains.
 * Created by bittrich on 7/10/17.
 */
public class UniProtFeatureContainer extends FeatureContainerEntry {
    private final Map<String, Feature> features;

    public UniProtFeatureContainer(AbstractFeatureProvider featureProvider) {
        super(featureProvider);
        this.features = new HashMap<>();
    }

    public void addFeature(String accession, Feature feature) {
        features.put(accession, feature);
    }

    public Map<String, Feature> getFeatures() {
        return features;
    }

    @Override
    public String toString() {
        return features.entrySet()
                .stream()
                .map(entry -> "@" + entry.getKey() + ": " + entry.getValue())
                .collect(Collectors.joining(", ", "[", "]"));
    }
}
