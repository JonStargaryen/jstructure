package design.aggregation;

import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.nio.file.Paths;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Move extracted fragments into a number of flexible clusters. Fragments are assigned to clusters when their respective
 * RMSD is <1 A, otherwise a new cluster is created. Clusters are significant, when they contain at least 20 entries.
 *
 * Diary:
 * <pre>
 *     12.14.16
 *      - rerun
 *      - reversed order: now most frequent clusters feature the lowest indices
 * </pre>
 * TODO extract sequences, create summary files, create folder hierarchy
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
                .forEach(motif -> handleMotif(topology, motif));
    }

    private static void handleMotif(String topology, String motif) {
        System.out.println("motif: " + motif);
        System.out.println("aligning structures");
        List<AtomContainer> proteins = DesignConstants.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                topology + "/"))
                .parallel()
                .filter(path -> path.getFileName().toString().startsWith(motif))
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());
        System.out.println("loaded " + proteins.size() + " atom containers");

        // compose clusters
        FragmentClusteringComposer fragmentClusteringComposer = new FragmentClusteringComposer(DesignConstants.IDENTICAL_CLUSTER_RMSD_THRESHOLD);
        fragmentClusteringComposer.composeClusterRepresentation(proteins);
        List<StructureCluster> clusters = fragmentClusteringComposer.getClusters();

        // sort clusters by their size
        System.out.println("composed " + clusters.size() + " clusters");
        clusters.sort(Comparator.comparingInt(o -> ((StructureCluster) o).getOriginalEntries().size()).reversed());
        System.out.println("cluster sizes: " + clusters.stream()
                .map(StructureCluster::getOriginalEntries)
                .mapToInt(Collection::size)
                .mapToObj(String::valueOf)
                .collect(Collectors.joining(", ", "[", "]")));

        // merge sparsely populated clusters
        List<StructureCluster> rareClusters = clusters.stream()
                .filter(cluster -> (double) cluster.getOriginalEntries().size() / (double) proteins.size() < DesignConstants.RARE_CLUSTER_THRESHOLD)
                .collect(Collectors.toList());
        System.out.println("rare clusters: " + rareClusters);
        if(rareClusters.size() != 0) {
            clusters.removeAll(rareClusters);
            //TODO really 'brute'-merge all lone-wolf-clusters into one consensus?
            StructureCluster mergedCluster = rareClusters.remove(0);
            rareClusters.stream()
                    .map(StructureCluster::getOriginalEntries)
                    .flatMap(Collection::stream)
                    .forEach(mergedCluster::add);
            clusters.add(0, mergedCluster);
        } else {
            clusters.add(0, StructureCluster.EMPTY_CLUSTER);
        }

        System.out.println("cluster sizes: " + clusters.stream()
                .map(StructureCluster::getOriginalEntries)
                .map(Collection::size)
                .map(String::valueOf)
                .collect(Collectors.joining(", ", "[", "]")));

        // output consensus
        for (int i = 0; i < clusters.size(); i++) {
            if(i == 0 && clusters.get(0).getOriginalEntries().size() == 0) {
                // no rare clusters were merged - thus, no freak cluster at index 0
                continue;
            }

            StructureCluster structureCluster = clusters.get(i);
            AtomContainer consensus = structureCluster.getConsensusRepresentation();
            DesignConstants.write(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR +
                    topology + "/" + motif + "-" + i + "-" + structureCluster.getOriginalEntries().size()
                    + DesignConstants.PDB_SUFFIX), consensus.composePDBRecord().getBytes());
        }
    }
}
