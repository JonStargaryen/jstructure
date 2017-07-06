package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the eGOR implementation.
 * Created by bittrich on 1/23/17.
 */
public class EnergyProfilePredictorTest {
    private Protein protein;
    private AbstractFeatureProvider energyProfilePredictor;

    @Before
    public void setup() {
        protein = ProteinParser.source("1brr").parse();
        energyProfilePredictor = FeatureProviderRegistry.resolvePredictor(EnergyProfile.class);
    }

    @Test
    public void shouldPredictEnergyProfile() {
        energyProfilePredictor.process(protein);

        protein.aminoAcids()
                .map(Group::getFeatureContainer)
                .map(container -> container.getFeature(EnergyProfile.class))
                .map(EnergyProfile::getSolvationEnergy)
                .forEach(value -> Assert.assertTrue(!value.isNaN()));
    }
}
