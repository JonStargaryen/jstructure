package alignment;

import de.bioforscher.jstructure.alignment.svd.ConsensusTreeComposer;
import de.bioforscher.jstructure.model.BinaryTree;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import static util.TestUtils.FRAGMENT_DIR;

/**
 * Ensures the {@link ConsensusTreeComposer} creates correct trees.
 * Created by S on 15.11.2016.
 */
public class ConsensusTreeComposerFunctionalTest {
    @Test
    public void shouldBuildConsensusTree() throws IOException {
        List<Protein> alignedFragments = Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
                .map(ProteinParser::parsePDBFile)
                .collect(StructureCollectors.toAlignedEnsemble());

        ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
        consensusTreeComposer.composeConsensusTree(alignedFragments);
        BinaryTree consensusTree = consensusTreeComposer.getConsensusTree();

        System.out.println(consensusTree.composeNewickRepresentation());
    }
}
