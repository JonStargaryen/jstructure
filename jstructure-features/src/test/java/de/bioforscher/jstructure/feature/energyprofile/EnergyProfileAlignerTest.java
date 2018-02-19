package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks the functionality of the {@link EnergyProfileAligner}.
 * Created by bittrich on 1/23/17.
 */
public class EnergyProfileAlignerTest {
    private Structure protein1acj;
    private Structure protein1brr;
    private EnergyProfileAligner energyProfileAligner;

    @Before
    public void setup() {
        protein1acj = StructureParser.fromPdbId("1acj").parse();
        protein1brr = StructureParser.fromPdbId("1brr").parse();
        FeatureProvider energyProfileCalculator = new EnergyProfileCalculator();
        energyProfileCalculator.process(protein1acj);
        energyProfileCalculator.process(protein1brr);
        energyProfileAligner = new EnergyProfileAligner();
    }

    @Test
    public void shouldResultInOptimalAlignment() {
        double distanceScore = energyProfileAligner.align(protein1acj, protein1acj);
        Assert.assertEquals(0.0, distanceScore, 0.0);
    }

    @Test
    public void shouldResultInWorstAlignment() {
        double distanceScore = energyProfileAligner.align(protein1acj, protein1brr);
        Assert.assertTrue(distanceScore > 0.0 && distanceScore < 5.0);
    }
}
