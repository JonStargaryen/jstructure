package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Group;

/**
 * Energy values associated to {@link Group} instances within a protein entity.
 * Created by bittrich on 5/17/17.
 */
@DefaultFeatureProvider(EnergyProfileCalculator.class)
public class EnergyProfile extends FeatureContainerEntry {
    private final double solvationEnergy;

    EnergyProfile(FeatureProvider featureProvider, double solvationEnergy) {
        super(featureProvider);
        this.solvationEnergy = solvationEnergy;
    }

    public double getSolvationEnergy() {
        return solvationEnergy;
    }
}
