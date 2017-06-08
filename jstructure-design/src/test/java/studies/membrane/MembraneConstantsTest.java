package studies.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.SecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Test capabilities of the MembraneConstants class.
 * Created by bittrich on 6/7/17.
 */
public class MembraneConstantsTest {
    @Test
    @Ignore("long creation of PDBTM data set skipped")
    public void getChainsOfPdbtmAlphaNrList() throws Exception {
        MembraneConstants.PdbtmAlphaNr.getChains()
                .forEach(chain -> {
                    System.out.println(chain.getChainId().getFullName());
                    chain.getFeatureContainer().getFeature(PLIPInteractionContainer.class);
                    chain.getFeatureContainer().getFeature(MembraneContainer.class);
                    chain.aminoAcids().findFirst().get().getFeatureContainer().getFeature(SecondaryStructure.class);
                });
    }
}