package design.aggregation;

import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * The table contains trans topology which is not needed and lacks information on the fragment size which eases
 * visualization.
 * Created by S on 23.11.2016.
 */
@Deprecated
public class S09_CleanUpFragmentStatistics {
    public static void main(String[] args) throws IOException {
        String output = Files.lines(Paths.get(DesignConstants.STATISTICS_DIR + "aligned-fragment-statistics.csv"))
                .filter(line -> !line.contains("trans"))
                .map(S09_CleanUpFragmentStatistics::handleLine)
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator()));

        Files.write(Paths.get(DesignConstants.STATISTICS_DIR + "aligned-fragment-statistics-reduced.csv"), output.getBytes());
    }

    private static String handleLine(String line) {
        if(line.startsWith("id")) {
            return line + ",size";
        } else {
            return line + "," + line.split("-")[0].substring(2);
        }
    }
}
