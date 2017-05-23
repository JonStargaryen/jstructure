package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Group;

/**
 * Energy values associated to {@link Group} instances within a protein entity.
 * Created by bittrich on 5/17/17.
 */
public class EnergyProfile extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double solvationEnergy;

    public EnergyProfile(AbstractFeatureProvider featureProvider, double solvationEnergy) {
        super(featureProvider);
        this.solvationEnergy = solvationEnergy;
    }

    public double getSolvationEnergy() {
        return solvationEnergy;
    }

    @Override
    public Double getValue() {
        return getSolvationEnergy();
    }
}
