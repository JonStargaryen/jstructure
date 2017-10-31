package de.bioforscher.jstructure.feature.evolution;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Map;

@DefaultFeatureProvider(EvolutionaryInformationCalculator.class)
public class EvolutionaryInformation extends FeatureContainerEntry {
    private final Map<AminoAcid.Family, Double> exchangeScores;
    private final double information;

    public EvolutionaryInformation(FeatureProvider featureProvider,
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

    public double getEvolutionaryPenality(AminoAcid.Family original, AminoAcid.Family mutation) {
        return exchangeScores.get(original) - exchangeScores.get(mutation);
    }
}
