package alignment;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static util.TestUtils.FRAGMENT_DIR;

/**
 * Created by S on 06.12.2016.
 */
public class FragmentClusteringComposerTest {
    @Test
    public void shouldBuildClusters() throws IOException {
        List<Protein> fragments = Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());

        FragmentClusteringComposer fragmentClusteringComposer = new FragmentClusteringComposer();
        fragmentClusteringComposer.composeClusterRepresentation(fragments);

        // will result in 1 cluster
        Assert.assertTrue(fragmentClusteringComposer.getClusters().size() == 2);
        AtomContainer consensus = fragmentClusteringComposer.getClusters().get(0).getConsensusRepresentation();

        SVDSuperimposer svdSuperimposer = new SVDSuperimposer();

        // which should be similar to each initial fragment
        double maxRmsd = fragments.stream()
                .map(fragment -> svdSuperimposer.align(fragment, consensus))
                .mapToDouble(AlignmentResult::getRmsd)
                .peek(rmsd -> System.out.println("rmsd of fragment to consensus: " + rmsd))
                .max()
                .orElseThrow(IllegalArgumentException::new);

        System.out.println("max rmsd: " + maxRmsd);
    }
}