package design.learning;

import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import static design.DesignConstants.DELIMITER;

/**
 * Depict clusters by WebLogos.
 * Created by bittrich on 12/19/16.
 */
public class S02_ConvertToWebLogo {
    private static final String motif = "GG4";
    private static final String topology = "ntm";
    private static final String clusterid = "4";

    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(DesignConstants.FRAGMENT_CLUSTERS_DIR + topology + "/" + motif + "/" + DesignConstants.CLUSTER_SUMMARY))
                .filter(line -> line.startsWith(clusterid))
                .map(line -> line.split(DELIMITER))
                .map(split -> split[4])
                .forEach(System.out::println);
    }
}
