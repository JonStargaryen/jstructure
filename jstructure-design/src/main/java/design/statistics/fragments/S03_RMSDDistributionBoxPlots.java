package design.statistics.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringJoiner;

import static design.DesignConstants.DELIMITER;
import static design.statistics.fragments.S02_PieChartsForTopologies.mapToTopology;

/**
 * Create box plots for the distribution of RMSD values for various models.
 * Created by bittrich on 12/20/16.
 */
public class S03_RMSDDistributionBoxPlots {
    public static void main(String[] args) throws IOException {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator(),
                "id" + DELIMITER + "motif" + DELIMITER + "topology" + DELIMITER + "rmsd" + DELIMITER + "size" + System.lineSeparator(),
                "");
        for (String topology : DesignConstants.TOPOLOGIES) {
            String mappedTopology = mapToTopology(topology);
            for (SequenceMotifDefinition sequenceMotifDefinition : SequenceMotifDefinition.values()) {
                Path summaryPath = Paths.get(DesignConstants.FRAGMENT_CLUSTERS_DIR + topology + "/" +
                        sequenceMotifDefinition.name() + "/" + DesignConstants.CLUSTER_SUMMARY);

                int size = Integer.valueOf(sequenceMotifDefinition.name().substring(2));

                Files.lines(summaryPath)
                        .filter(line -> !line.startsWith("clusterId"))
                        .map(line -> line.split(DELIMITER))
                        .map(split -> split[3] + DELIMITER + sequenceMotifDefinition.name() + DELIMITER + mappedTopology +
                                DELIMITER + split[5] + DELIMITER + size)
                        .peek(System.out::println)
                        .forEach(stringJoiner::add);
            }
        }

        Files.write(Paths.get(DesignConstants.TSV_DIR + "rmsd-distribution.tsv"), stringJoiner.toString().getBytes());
    }
}
