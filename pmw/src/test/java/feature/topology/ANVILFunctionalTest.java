package feature.topology;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests the ANVIL algorithm which screens for the most suitable membrane orientation for a given protein.
 * Created by S on 03.11.2016.
 */
public class ANVILFunctionalTest {
    private Protein protein1brr;

    @Before
    public void setup() {
        //TODO algorithms should compute needed features automatically
        protein1brr = ProteinParser.parseProteinById("1brr");
        new AccessibleSurfaceAreaCalculator().process(protein1brr);
    }

    @Test
    public void shouldRunANVIL() {
        //TODO test
        new ANVIL().process(protein1brr);
        System.out.println("membrane quality is " + protein1brr.getFeature(Membrane.class,
                ANVIL.FeatureNames.MEMBRANE).getQmax());
    }
}