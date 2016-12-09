package design.statistics;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static design.DesignConstants.DELIMITER;

/**
 * Count occurrence of each fragment cluster. This time however, fragments occurring less than 5 times are merged into a
 * 'others' cluster.
 * Created by S on 08.12.2016.
 */
public class S06_PieChartOccurrenceReduced {
    private static final String TOPOLOGY = "trans";
    private static final int MIN_OCCURRENCE = 20;

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

                    List<String> lines = DesignConstants.list(topologyDir)
                            .filter(path -> path.toFile().getName().startsWith(motif.name()))
                            // skip entries of low frequency
                            .filter(path -> Integer.parseInt(path.toFile().getName().split("-")[2].split("\\.")[0]) >= MIN_OCCURRENCE)
                            .map(path -> motif.name() + DELIMITER + count + DELIMITER +
                                    path.toFile().getName().split("-")[1] + DELIMITER +
                                    path.toFile().getName().split("-")[2].split("\\.")[0])
                            .collect(Collectors.toList());

                    int occurenceOfRareClusters = DesignConstants.list(topologyDir)
                            .filter(path -> path.toFile().getName().startsWith(motif.name()))
                            // skip entries of low frequency
                            .filter(path -> Integer.parseInt(path.toFile().getName().split("-")[2].split("\\.")[0]) < MIN_OCCURRENCE)
                            .map(path -> path.toFile().getName().split("-")[2].split("\\.")[0])
                            .mapToInt(Integer::valueOf)
                            .sum();

                    lines.add(motif.name() + DELIMITER + count + DELIMITER + "0" + DELIMITER + occurenceOfRareClusters);

                    // renumber clusters
                    List<String> renumberedLines = new ArrayList<>();
                    int clusterNumber = 1;
                    String clusterMotif = "";
                    for(String line : lines) {
                        String[] tmp = line.split(DELIMITER);
                        if(tmp[2].equals("0")) {
                            renumberedLines.add(line);
                            continue;
                        }

                        // motif changed
                        if(!tmp[0].equals(clusterMotif)) {
                            clusterMotif = tmp[0];
                            clusterNumber = 1;
                        }

                        renumberedLines.add(tmp[0] + DELIMITER + tmp[1] + DELIMITER + clusterNumber + DELIMITER + tmp[3]);
                        clusterNumber++;
                    }

                    return renumberedLines.stream();
                })
                .collect(Collectors.joining(System.lineSeparator(), "motif" + DELIMITER + "abs" + DELIMITER +
                        "cluster" + DELIMITER + "count" + System.lineSeparator(), ""));
    }
}

