package alignment;

import de.bioforscher.jstructure.alignment.consensus.ConsensusTreeComposer;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;
import util.TestUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static util.TestUtils.FRAGMENT_DIR;

/**
 * Ensures the {@link ConsensusTreeComposer} creates correct trees.
 * Created by S on 15.11.2016.
 */
public class ConsensusTreeComposerFunctionalTest {
    @Test
    public void shouldBuildConsensusTree() throws IOException {
        List<Protein> fragments = Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());

        ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
        consensusTreeComposer.composeConsensusTree(fragments);

        System.out.println(consensusTreeComposer.getConsensusTree().composeNewickRepresentation());
        consensusTreeComposer.getAlignedContainers().stream()
                .map(AtomContainer::composePDBRecord)
                .map(string -> string + System.lineSeparator())
                .forEach(System.out::println);
    }

    @Test
    public void shouldMergeTwoFragments() throws IOException {
        List<Protein> fragments = Files.list(Paths.get(TestUtils.getResourceAsFilepath(FRAGMENT_DIR)))
                // skip 10 entries, last two are most dissimilar
                .skip(10)
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());

        ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
        consensusTreeComposer.composeConsensusTree(fragments);

        System.out.println(consensusTreeComposer.getConsensusTree().composeNewickRepresentation());
        consensusTreeComposer.getAlignedContainers().stream()
                .map(AtomContainer::composePDBRecord)
                .map(string -> string + System.lineSeparator())
                .forEach(System.out::println);
    }

    @Test
    public void shouldSerializeConsensusTree() {

    }
}
