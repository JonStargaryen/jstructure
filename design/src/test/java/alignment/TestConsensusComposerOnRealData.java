package alignment;

import de.bioforscher.jstructure.alignment.consensus.ConsensusTreeComposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Tests the {@link ConsensusTreeComposer} on real-world data.
 * Created by S on 23.11.2016.
 */
public class TestConsensusComposerOnRealData {
    private SequenceMotifDefinition motif;

    @Before
    public void setup() {
        this.motif = SequenceMotifDefinition.GG4;
    }

    @Test
    public void shouldComposeConsensusMotifForTMRegion() throws IOException {
        List<AtomContainer> alignedEnsemble = getAlignedEnsemble("tm");

        AtomContainer reference = alignedEnsemble.remove(0);
        double averageReferenceRmsd = alignedEnsemble.stream()
                .filter(protein -> reference.getAtoms().size() == protein.getAtoms().size())
                .mapToDouble(protein -> LinearAlgebraAtom.calculateRmsd(reference, protein))
                .average()
                .orElse(-1);

        ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
        consensusTreeComposer.composeConsensusTree(alignedEnsemble);
        AtomContainer consensus = consensusTreeComposer.getConsensus();

        double averageConsensusRmsd = alignedEnsemble.stream()
                .filter(protein -> consensus.getAtoms().size() == protein.getAtoms().size())
                .mapToDouble(protein -> LinearAlgebraAtom.calculateRmsd(consensus, protein))
                .average()
                .orElse(-1);

        System.out.println("reference size: " + reference.getAtoms().size());
        System.out.println("consensus size: " + consensus.getAtoms().size());
        System.out.println(reference.atoms()
                .map(Atom::getName)
                .collect(Collectors.joining()));
        System.out.println(consensus.atoms()
                .map(Atom::getName)
                .collect(Collectors.joining()));
        System.out.println("average reference RMSD for tm: " + averageReferenceRmsd);
        System.out.println("average consensus RMSD for tm: " + averageConsensusRmsd);
}

    @Test
    public void shouldComposeConsensusMotifForNTMRegion() throws IOException {
        List<AtomContainer> alignedEnsemble = getAlignedEnsemble("ntm");

        AtomContainer reference = alignedEnsemble.remove(0);
        double averageRmsd = alignedEnsemble.stream()
                .filter(protein -> reference.getAtoms().size() == protein.getAtoms().size())
                .limit(100)
                .mapToDouble(protein -> LinearAlgebraAtom.calculateRmsd(reference, protein))
                .average()
                .getAsDouble();

        System.out.println("average RMSD for ntm: " + averageRmsd);
    }

    private List<AtomContainer> getAlignedEnsemble(String topology) throws IOException {
        return Files.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/"))
                .filter(path -> path.toFile().getName().startsWith(motif.name()))
                .limit(100)
                .map(ProteinParser::parsePDBFile)
                .collect(StructureCollectors.toAlignedEnsembleByConsensus());
    }
}
