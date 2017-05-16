package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
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
        energyProfilePredictor = FeatureProviderRegistry.resolvePredictor(EnergyProfileCalculator.SOLVATION_ENERGY);
    }

    @Test
    public void shouldPredictEnergyProfile() {
        energyProfilePredictor.process(protein);

        protein.aminoAcids()
                .map(group -> group.getFeatureAsDouble(EnergyProfileCalculator.SOLVATION_ENERGY))
                .forEach(value -> Assert.assertTrue(!value.isNaN()));
    }
}
