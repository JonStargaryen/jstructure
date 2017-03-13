package design.statistics;

import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import static design.DesignConstants.DELIMITER;

/**
 * Compute statistics on sets of aligned fragments grouped by topology.
 * Created by S on 11.11.2016.
 */
@Deprecated
public class S01_AlignedFragmentStatistics {
    public static void main(String[] args) throws IOException {
        StringBuilder output = new StringBuilder();
        output.append("id")
                .append(DELIMITER)
                .append("motif")
                .append(DELIMITER)
                .append("topo")
                .append(DELIMITER)
                .append("clusterIndex")
                .append(DELIMITER)
                .append("rmsd")
                .append(System.lineSeparator());

        for(SequenceMotifDefinition definition : SequenceMotifDefinition.values()) {
            Files.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR))
                    .forEach((Path topologyDir) -> {
                            System.out.println(topologyDir);
                            // align ensemble
//                                List<AtomContainer> alignedEnsemble = Files.list(topologyDir)
//                                        .filter(path -> path.toFile().getName().startsWith(definition.name()))
//                                        .map(path -> ProteinParser.source(path).parse())
//                                        .collect(StructureCollectors.toAlignedEnsembleByConsensus());
//
//                                // create consensus
//                                ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
//                                consensusTreeComposer.composeConsensusTree(alignedEnsemble);
//                                AtomContainer consensus = consensusTreeComposer.getConsensus();

//                        for(AtomContainer observation : alignedEnsemble) {
//                            double rmsd = LinearAlgebraAtom.calculateRmsd(Selection.on(observation)
//                                            .backboneAtoms()
//                                            .asAtomContainer(),
//                                    Selection.on(consensus)
//                                            .backboneAtoms()
//                                            .asAtomContainer());
//                            String line = observation.getIdentifier() + delimiter + definition.name() +
//                                    delimiter + topologyDir.toFile().getName() + delimiter + rmsd + System.lineSeparator();
//                            System.out.println(line);
//                            output.append(line);
//                        }

                            // align ensemble
                            List<StructureCluster> clusters = DesignConstants.list(topologyDir)
                                    .filter(path -> path.toFile().getName().startsWith(definition.name()))
                                    .map(path -> ProteinParser.source(path).parse())
                                    .collect(StructureCollectors.toAlignedEnsembleByStructuralClustering());

                            for(int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++) {
                                StructureCluster cluster = clusters.get(clusterIndex);
                                AtomContainer consensus = cluster.getConsensusRepresentation();
                                for(AtomContainer entry : cluster.getOriginalEntries()) {
                                    double rmsd = LinearAlgebraAtom.calculateRmsd(consensus, cluster.getConsensusRepresentation());
                                    String line = entry.getIdentifier() + DELIMITER + definition.name() + DELIMITER +
                                            topologyDir.toFile().getName() + DELIMITER + clusterIndex + DELIMITER + rmsd
                                            + System.lineSeparator();
                                    System.out.print(line);
                                    output.append(line);
                                }
                            }
                    });

            Files.write(Paths.get(DesignConstants.STATISTICS_DIR + "aligned-fragment-statistics.csv"), output.toString().getBytes());
        }
    }
}
