package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

@DefaultFeatureProvider(EgorAgreementCalculator.class)
public class EgorAgreement extends EnergyProfile {
    private final double egorPrediction;
    private final double egorAgreement;

    public EgorAgreement(FeatureProvider featureProvider,
                         double solvationEnergy,
                         double egorPrediction,
                         double egorAgreement) {
        super(featureProvider, solvationEnergy);
        this.egorPrediction = egorPrediction;
        this.egorAgreement = egorAgreement;
    }

    public double getEgorPrediction() {
        return egorPrediction;
    }

    public double getEgorAgreement() {
        return egorAgreement;
    }
}
