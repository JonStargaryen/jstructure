package design.aggregation;

import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Move extracted fragments into a number of clusters
 * Created by S on 06.12.2016.
 */
public class S10_BuildFragmentClusters {
    public static void main(String[] args) {
        S06_ExtractSequences.TOPOLOGIES.stream()
                .skip(1)
                .forEach(S10_BuildFragmentClusters::handleTopology);
    }

    private static void handleTopology(final String topology) {
        System.out.println("topology: " + topology);
        Stream.of(SequenceMotifDefinition.values())
                .skip(19)
                .map(SequenceMotifDefinition::name)
                .forEach((String motif) -> {
                    System.out.println("motif: " + motif);
                    try {
                        System.out.println("aligning structures");
                        List<AtomContainer> proteins = Files.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                                topology + "/"))
                                .parallel()
                                .filter(path -> path.getFileName().toString().startsWith(motif))
                                .map(ProteinParser::parsePDBFile)
                                .collect(Collectors.toList());

                        FragmentClusteringComposer fragmentClusteringComposer = new FragmentClusteringComposer();
                        fragmentClusteringComposer.composeClusterRepresentation(proteins);

                        List<FragmentClusteringComposer.StructureCluster> clusters = fragmentClusteringComposer.getClusters();
                        System.out.println("composed " + clusters.size() + " clusters");

                        for (int i = 0; i < clusters.size(); i++) {
                            FragmentClusteringComposer.StructureCluster structureCluster = clusters.get(i);
                            AtomContainer consensus = structureCluster.getConsensusRepresentation();
                            Files.write(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR +
                                    topology + "/" + motif + "-" + (i + 1) + "-" + structureCluster.getEntries().size() + DesignConstants.PDB_SUFFIX), consensus.composePDBRecord().getBytes());
                        }
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}
