package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the eGOR implementation.
 * Created by bittrich on 1/23/17.
 */
public class EnergyProfilePredictorTest {
    private Structure protein;
    private FeatureProvider energyProfilePredictor;

    @Before
    public void setup() {
        protein = StructureParser.source("1brr").parse();
        energyProfilePredictor = new EnergyProfilePredictor();
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
