package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Information obtained by a UniProt mapping to homologous protein chains.
 * Created by bittrich on 7/10/17.
 */
@DefaultFeatureProvider(UniProtHomologyAnnotator.class)
public class UniProtFeatureContainer extends FeatureContainerEntry {
    private final Map<String, List<Feature>> features;

    public UniProtFeatureContainer(FeatureProvider featureProvider) {
        super(featureProvider);
        this.features = new HashMap<>();
    }

    public void addFeature(String accession, Feature feature) {
        if(!features.containsKey(accession)) {
            features.put(accession, new ArrayList<>());
        }
        List<Feature> featureList = features.get(accession);
        featureList.add(feature);
    }

    public Map<String, List<Feature>> getFeatureMap() {
        return features;
    }

    public List<Feature> getFeatures(FeatureType... featureTypes) {
        List<FeatureType> featureTypeList = Arrays.asList(featureTypes);
        return features.values()
                .stream()
                .flatMap(Collection::stream)
                .filter(feature -> featureTypeList.contains(feature.getType()))
                .collect(Collectors.toList());
    }

    @Override
    public String toString() {
        return features.entrySet()
                .stream()
                .map(entry -> "@" + entry.getKey() + ": " + entry.getValue())
                .collect(Collectors.joining(", ", "[", "]"));
    }
}
