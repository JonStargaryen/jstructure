package de.bioforscher.jstructure.feature.evolution;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Map;

public class EvolutionaryInformation extends FeatureContainerEntry {
    private final Map<AminoAcid.Family, Double> exchangeScores;
    private final double information;

    public EvolutionaryInformation(AbstractFeatureProvider featureProvider,
                                   Map<AminoAcid.Family, Double> exchangeScores,
                                   double information) {
        super(featureProvider);
        this.exchangeScores = exchangeScores;
        this.information = information;
    }

    public Map<AminoAcid.Family, Double> getExchangeScores() {
        return exchangeScores;
    }

    public double getInformation() {
        return information;
    }
}
