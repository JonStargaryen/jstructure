package design.aggregation;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Stream;

/**
 * Aligns fragments of one sequence motif observation.
 * Deprecated because the consensus-by-reference approach was replaced by the {@link FragmentClusteringComposer}.
 * Created by S on 07.11.2016.
 */
@Deprecated
public class S07_AlignFragments {
    public static void main(String[] args) {
        S06_ExtractSequences.TOPOLOGIES.forEach(S07_AlignFragments::handleTopology);
    }

    private static void handleTopology(final String topology) {
        System.out.println("topology: " + topology);
        Stream.of(SequenceMotifDefinition.values())
                .map(SequenceMotifDefinition::name)
                .forEach(motif -> {
                        System.out.println("motif: " + motif);
                        try {
                            System.out.println("aligning structures");
                            List<AtomContainer> alignedProteins = Files.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                                            topology + "/"))
                                    .filter(path -> path.getFileName().toString().startsWith(motif))
                                    .map(ProteinParser::parsePDBFile)
                                    .collect(StructureCollectors.toAlignedEnsembleByConsensus());

                            System.out.println("writing aligned files");
                            // write structures
                            alignedProteins.forEach(protein -> {
                                try {
                                    Files.write(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                                            topology + "/" + protein.getIdentifier() + "-" + protein.getFeatureAsDouble(AbstractAlignmentAlgorithm.FRAGMENT_RMSD)
                                                    + DesignConstants.PDB_SUFFIX),
                                            protein.composePDBRecord().getBytes());
                                } catch (IOException e) {
                                    throw new UncheckedIOException(e);
                                }
                            });
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    });
    }
}
