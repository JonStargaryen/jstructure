package design;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;
import java.util.Map;

/**
 * The test cases for the ProteinSource.
 * Created by S on 13.12.2016.
 */
public class ProteinSourceFunctionalTest {
    @Test
    public void shouldRetrieveAllFragmentClusters() {
        Map<String, List<AtomContainer>> fragments = ProteinSource.getGroupedFragments("tm", SequenceMotifDefinition.LY6);

        Assert.assertNotNull(fragments);

        System.out.println("number of clusters: " + fragments.size());

        fragments.entrySet().forEach(System.out::println);
    }
}
