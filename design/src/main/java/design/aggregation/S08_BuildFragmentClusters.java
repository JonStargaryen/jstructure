package design.aggregation;

import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Move extracted fragments into a number of flexible clusters. Fragments are assigned to clusters when their respective
 * RMSD is <1 A, otherwise a new cluster is created. Clusters are significant, when they contain at least 20 entries.
 * Created by S on 06.12.2016.
 */
public class S08_BuildFragmentClusters {
    public static void main(String[] args) {
        DesignConstants.TOPOLOGIES.forEach(S08_BuildFragmentClusters::handleTopology);
    }

    private static void handleTopology(final String topology) {
        System.out.println("topology: " + topology);
        Stream.of(SequenceMotifDefinition.values())
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
                        System.out.println("loaded " + proteins.size() + " atom containers");

                        FragmentClusteringComposer fragmentClusteringComposer = new FragmentClusteringComposer();
                        fragmentClusteringComposer.composeClusterRepresentation(proteins);

                        List<StructureCluster> clusters = fragmentClusteringComposer.getClusters();
                        System.out.println("composed " + clusters.size() + " clusters");
                        clusters.sort(Comparator.comparingInt(o -> o.getEntries().size()));
                        System.out.println("cluster sizes: " + clusters.stream()
                                .map(StructureCluster::getEntries)
                                .mapToInt(Collection::size)
                                .mapToObj(String::valueOf)
                                .collect(Collectors.joining(", ", "[", "]")));

                        for (int i = 0; i < clusters.size(); i++) {
                            StructureCluster structureCluster = clusters.get(i);
                            AtomContainer consensus = structureCluster.getConsensusRepresentation();
                            Files.write(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR +
                                    topology + "/" + motif + "-" + (i + 1) + "-" + structureCluster.getEntries().size()
                                    + DesignConstants.PDB_SUFFIX), consensus.composePDBRecord().getBytes());
                        }
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}
