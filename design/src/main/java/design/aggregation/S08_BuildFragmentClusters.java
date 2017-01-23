package design.aggregation;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.nio.file.Paths;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static design.DesignConstants.DELIMITER;

/**
 * Move extracted fragments into a number of flexible clusters. Fragments are assigned to clusters when their respective
 * RMSD is <1 A, otherwise a new cluster is created. Clusters are significant, when they contain at least 20 entries.
 *
 * Diary:
 * <pre>
 *     12/14/16
 *      - reversed order: now most frequent clusters feature the lowest indices
 *     12/15/16
 *      - writes summary file
 *      - writes aligned fragments
 *      - rerun with cutoff < 1 A
 *          - this could be too generous -> rerun with more rigorous threshold for tm/trans topology
 * </pre>
 * Created by S on 06.12.2016.
 */
public class S08_BuildFragmentClusters {
    public static void main(String[] args) {
        DesignConstants.TOPOLOGIES.stream()
                // skip ntm for now
                .filter(topology -> !topology.equals("ntm"))
                .forEach(S08_BuildFragmentClusters::handleTopology);
    }

    private static void handleTopology(final String topology) {
        System.out.println("topology: " + topology);
        Stream.of(SequenceMotifDefinition.values())
                .map(SequenceMotifDefinition::name)
                .forEach(motif -> handleMotif(topology, motif));
    }

    private static void handleMotif(String topology, String motif) {
        String basePath = DesignConstants.FRAGMENT_CLUSTERS_DIR + topology + "/" + motif + "/";
        DesignConstants.makeDirectoryIfAbsent(Paths.get(basePath));

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
        FragmentClusteringComposer fragmentClusteringComposer = new FragmentClusteringComposer(DesignConstants.IDENTICAL_CLUSTER_RMSD_THRESHOLD, false);
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
        StringBuilder output = new StringBuilder();
        output.append("clusterId")
                .append(DELIMITER)
                .append("totalCount")
                .append(DELIMITER)
                .append("clusterSize")
                .append(DELIMITER)
                .append("fragmentId")
                .append(DELIMITER)
                .append("fragmentSequence")
                .append(DELIMITER)
                .append("rmsd")
                .append(System.lineSeparator());

        for (int i = 0; i < clusters.size(); i++) {
            if(i == 0 && clusters.get(0).getOriginalEntries().size() == 0) {
                // no rare clusters were merged - thus, no freak cluster at index 0
                continue;
            }

            String clusterPath = basePath + i + "/";
            DesignConstants.makeDirectoryIfAbsent(Paths.get(clusterPath));

            StructureCluster structureCluster = clusters.get(i);
            AtomContainer consensus = structureCluster.getConsensusRepresentation();
            DesignConstants.write(Paths.get(clusterPath + DesignConstants.CLUSTER_CONSENSUS),
                    consensus.composePDBRecord().getBytes());

            for(AtomContainer fragment : structureCluster.getOriginalEntries()) {
                // align fragment relative to consensus
                SVDSuperimposer svdSuperimposer = new SVDSuperimposer();
                AlignmentResult alignmentResult = svdSuperimposer.align(consensus, fragment);
                String rmsd = DesignConstants.DECIMAL_FORMAT.format(alignmentResult.getAlignmentScore());

                // add to summary file
                output.append(i)
                        .append(DELIMITER)
                        .append(proteins.size())
                        .append(DELIMITER)
                        .append(structureCluster.getOriginalEntries().size())
                        .append(DELIMITER)
                        .append(fragment.getIdentifier())
                        .append(DELIMITER)
                        .append(((GroupContainer) fragment).getAminoAcidSequence())
                        .append(DELIMITER)
                        .append(rmsd)
                        .append(System.lineSeparator());

                // write aligned file
                DesignConstants.write(Paths.get(clusterPath + fragment.getIdentifier() + DesignConstants.PDB_SUFFIX),
                        alignmentResult.getAlignedQuery().composePDBRecord().getBytes());
            }
        }

        DesignConstants.write(Paths.get(basePath + DesignConstants.CLUSTER_SUMMARY), output.toString().getBytes());
    }
}
