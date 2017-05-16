package design.statistics;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.stream.Collectors;

import static design.DesignConstants.DELIMITER;

/**
 * Count occurrence of each fragment cluster.
 * Created by S on 08.12.2016.
 */
@Deprecated
public class S05_PieChartOccurrence {
    private static final String TOPOLOGY = "tm";

    public static void main(String[] args) {
        String output = countOccurrences(TOPOLOGY);
        System.out.println(output);
    }

    private static String countOccurrences(String topology) {
        return Arrays.stream(SequenceMotifDefinition.values())
                .flatMap(motif -> {
                    Path topologyDir = Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR + topology + "/");

                    int count = DesignConstants.list(topologyDir)
                            .filter(path -> path.toFile().getName().startsWith(motif.name()))
                            .map(Path::toFile)
                            .map(File::getName)
                            .map(name -> name.split("-")[2].split("\\.")[0])
                            .mapToInt(Integer::valueOf)
                            .sum();

                    return DesignConstants.list(topologyDir)
                            .filter(path -> path.toFile().getName().startsWith(motif.name()))
                            .map(path -> motif.name() + DELIMITER + count + DELIMITER +
                                    path.toFile().getName().split("-")[1] + DELIMITER +
                                    path.toFile().getName().split("-")[2].split("\\.")[0]);
                })
                .collect(Collectors.joining(System.lineSeparator(), "motif" + DELIMITER + "abs" + DELIMITER + "cluster" + DELIMITER + "count" + System.lineSeparator(), ""));
    }
}
