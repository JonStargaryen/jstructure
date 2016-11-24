package design.aggregation;

import de.bioforscher.jstructure.alignment.svd.ConsensusTreeComposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.mathematics.CoordinateManipulations;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

/**
 * Compute statistics on sets of aligned fragments grouped by topology.
 * Created by S on 11.11.2016.
 */
public class S08_AlignedFragmentStatistics {
    public static void main(String[] args) throws IOException {
        StringBuilder output = new StringBuilder();
        final String delimiter = ",";
        output.append("id")
                .append(delimiter)
                .append("motif")
                .append(delimiter)
                .append("topo")
                .append(delimiter)
                .append("rmsd")
                .append(System.lineSeparator());

        for(SequenceMotifDefinition definition : SequenceMotifDefinition.values()) {
            Files.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMNET_BY_TOPOLOGY_DIR))
                    .forEach(topologyDir -> {
                            System.out.println(topologyDir);
                            try {
                                // align ensemble
                                List<Protein> alignedEnsemble = Files.list(topologyDir)
                                        .filter(path -> path.toFile().getName().startsWith(definition.name()))
                                        .map(ProteinParser::parsePDBFile)
                                        .collect(StructureCollectors.toAlignedEnsemble());

                                // create consensus
                                ConsensusTreeComposer consensusTreeComposer = new ConsensusTreeComposer();
                                consensusTreeComposer.composeConsensusTree(alignedEnsemble);
                                AtomContainer consensus = consensusTreeComposer.getConsensus();

                                for(AtomContainer observation : alignedEnsemble) {
                                    double rmsd = CoordinateManipulations.calculateRMSD(Selection.on(observation)
                                            .backboneAtoms()
                                            .asAtomContainer(),
                                            Selection.on(consensus)
                                            .backboneAtoms()
                                            .asAtomContainer());
                                    String line = observation.getIdentifier() + delimiter + definition.name() +
                                            delimiter + topologyDir.toFile().getName() + delimiter + rmsd + System.lineSeparator();
                                    System.out.println(line);
                                    output.append(line);
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                    });

            Files.write(Paths.get(DesignConstants.STATISTICS_DIR + "aligned-fragment-statistics.csv"), output.toString().getBytes());
        }
    }
}
