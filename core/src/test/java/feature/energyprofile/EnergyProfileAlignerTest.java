package feature.energyprofile;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileAligner;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Checks the functionality of the {@link EnergyProfileAligner}.
 * Created by bittrich on 1/23/17.
 */
public class EnergyProfileAlignerTest {
    private Protein protein1acj;
    private Protein protein1brr;
    private EnergyProfileAligner energyProfileAligner;

    @Before
    public void setup() {
        protein1acj = ProteinParser.parseProteinById("1acj");
        protein1brr = ProteinParser.parseProteinById("1brr");
        AbstractFeatureProvider energyProfileCalculator = FeatureProviderRegistry.resolve(EnergyProfileCalculator.SOLVATION_ENERGY);
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
