package design.statistics.fragments;

import design.DesignConstants;

import java.nio.file.Paths;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static design.DesignConstants.DELIMITER;

/**
 * Pie charts for the distribution of sequence motifs in various topologies.
 * Created by bittrich on 12/15/16.
 */
public class S02_PieChartsForTopologies {
    public static void main(String[] args) {
        System.out.println(Stream.of("tm", "ntm", "trans")
                .flatMap(topology -> {
                    String filepath = "/home/bittrich/" + topology + ".tsv";
                    return DesignConstants.lines(Paths.get(filepath))
                            .skip(1)
                            .map(line -> line.split(DELIMITER))
                            .map(split -> split[0] + DELIMITER + mapToTopology(topology) + DELIMITER + split[2]);
                })
                .distinct()
                .collect(Collectors.joining(System.lineSeparator())));
    }

    static String mapToTopology(String topology) {
        switch (topology) {
            case "tm":
                return "transmembrane";
            case "ntm":
                return "non-transmembrane";
            default:
                return "transition";
        }
    }
}
