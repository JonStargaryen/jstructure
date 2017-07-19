package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.*;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;

@FeatureProvider(provides = EnergyProfile.class)
public class EgorAgreementCalculator extends AbstractFeatureProvider {
    private final EnergyProfileCalculator energyProfileCalculator;
    private final EnergyProfilePredictor energyProfilePredictor;

    public EgorAgreementCalculator() {
        this.energyProfileCalculator = new EnergyProfileCalculator();
        this.energyProfilePredictor = new EnergyProfilePredictor();
    }

    @Override
    protected void processInternally(Structure protein) {
        FeatureContainer featureContainer = protein.select().aminoAcids().asAminoAcid().getFeatureContainer();
        // need to compute energy profile first
        if(featureContainer.getFeatures(EnergyProfile.class)
                .stream()
                .map(FeatureContainerEntry::getFeatureProvider)
                .noneMatch(energyProfileCalculator.getClass()::isInstance)) {
            energyProfileCalculator.process(protein);
        }
        // need to predict energy profile first
        if(featureContainer.getFeatures(EnergyProfile.class)
                .stream()
                .map(FeatureContainerEntry::getFeatureProvider)
                .noneMatch(energyProfilePredictor.getClass()::isInstance)) {
            energyProfilePredictor.process(protein);
        }

        protein.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        chain.aminoAcids().forEach(this::processInternally);
    }

    private void processInternally(AminoAcid aminoAcid) {
        List<EnergyProfile> energy = aminoAcid.getFeatureContainer().getFeatures(EnergyProfile.class);
        double calculation = energy.stream()
                .filter(energyProfile -> energyProfile.getFeatureProvider().getClass().isInstance(energyProfileCalculator))
                .findFirst()
                .map(EnergyProfile::getSolvationEnergy)
                .orElse(-1000.0);
        double prediction = energy.stream()
                .filter(energyProfile -> energyProfile.getFeatureProvider().getClass().isInstance(energyProfilePredictor))
                .findFirst()
                .map(EnergyProfile::getSolvationEnergy)
                .orElse(1000.0);
        double quantil = EnergyProfilePredictor.mapEnergyToQuantil(calculation);
        aminoAcid.getFeatureContainer().addFeature(new EgorAgreement(this,
                calculation,
                prediction,
                prediction == quantil ? 1.0 : 0.0));
    }
}
