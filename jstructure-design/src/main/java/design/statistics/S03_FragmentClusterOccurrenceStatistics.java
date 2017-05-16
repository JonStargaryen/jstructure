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
 * Starting point are TOPOLOGY > MOTIF > CLUSTERS.
 * Created by S on 06.12.2016.
 */
@Deprecated
public class S03_FragmentClusterOccurrenceStatistics {
    public static void main(String[] args) {
        System.out.println("motif" + DELIMITER + "topology" + DELIMITER + "count");
        Arrays.stream(SequenceMotifDefinition.values())
                .map(S03_FragmentClusterOccurrenceStatistics::countOccurrences)
                .forEach(System.out::println);
    }

    private static String countOccurrences(SequenceMotifDefinition sequenceMotifDefinition) {
        String motifName = sequenceMotifDefinition.name();
        return DesignConstants.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR))
                .map(topology -> {
                    int count = DesignConstants.list(topology)
                            .filter(path -> path.toFile().getName().startsWith(motifName))
                            .map(Path::toFile)
                            .map(File::getName)
                            .map(name -> name.split("-")[2].split("\\.")[0])
                            .mapToInt(Integer::valueOf)
                            .sum();
                    return motifName + DELIMITER + toTopologyString(topology) + DELIMITER + count;
                })
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String toTopologyString(Path topology) {
        String name = topology.toFile().getName();
        if(name.equals("tm")) {
            return "transmembrane";
        }
        if(name.equals("ntm")) {
            return "non-transmembrane";
        }
        return "transition";
    }
}
