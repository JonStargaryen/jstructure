package design.statistics.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static design.DesignConstants.DELIMITER;

/**
 * Convert the data to a TSV to generate the pie-chart-plot by R.
 * Created by bittrich on 12/15/16.
 */
public class S01_PieChartsForFragments {
    public static void main(String[] args) throws IOException {
        String topology = "ntm";

        String data = handleTopology(topology);
        System.out.println(data);
        Files.write(Paths.get("/home/bittrich/" + topology + ".tsv"), data.getBytes());
    }

    private static String handleTopology(String topology) {
        return Arrays.stream(SequenceMotifDefinition.values())
                .map(motif -> handleMotif(topology, motif))
                .collect(Collectors.joining(System.lineSeparator(),
                        "motif" + DELIMITER + "clusterId" + DELIMITER + "totalCount" + DELIMITER + "clusterSize" + System.lineSeparator(),
                        ""));
    }

    private static String handleMotif(String topology, SequenceMotifDefinition motif) {
        return DesignConstants.lines(Paths.get(DesignConstants.FRAGMENT_CLUSTERS_DIR + topology + "/" + motif + "/" + DesignConstants.CLUSTER_SUMMARY))
                //skip old header
                .filter(line -> !line.startsWith("clusterId"))
                .map(line -> motif + DELIMITER + Pattern.compile(DELIMITER).splitAsStream(line)
                                                                      .limit(3)
                                                                      .collect(Collectors.joining(DELIMITER))
                )
                .distinct()
                .collect(Collectors.joining(System.lineSeparator()));
    }
}
