package de.bioforscher.thermostability.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class A03_CreateLatexTemplate {
    private static String lastThermophilePdbId = "";
    private static String lastMesophilePdbId = "";
    private static List<String> featureNames;
    // values of one protein
    private static List<List<Double>> thermophileProtein = new ArrayList<>();
    private static List<List<Double>> mesophileProtein = new ArrayList<>();
    private static StringJoiner output = new StringJoiner(System.lineSeparator());

    public static void main(String[] args) throws IOException {
        List<String> lines = Files.readAllLines(Start2FoldConstants.THERMOSTABILITY_DIRECTORY
                .resolve("statistics")
                .resolve("thermophile.csv"));

        String[] headerLineSplit = lines.remove(0).split(",");
        featureNames = Arrays.stream(headerLineSplit)
                .skip(8)
                .limit(51)
                .collect(Collectors.toList());

        for(String line : lines) {
            String[] split = line.split(",");
            String pdbId = split[0];
            boolean thermophile = split[59].equals("thermophile");

            if(thermophile && !lastThermophilePdbId.equals(pdbId) && !lastThermophilePdbId.equals("")) {
                processPair();
            }

            List<Double> lineValues = Arrays.stream(split)
                    .skip(8)
                    .limit(51)
                    .map(Double::valueOf)
                    .collect(Collectors.toList());
            if(thermophile) {
                lastThermophilePdbId = pdbId;
                thermophileProtein.add(lineValues);
            } else {
                lastMesophilePdbId = pdbId;
                mesophileProtein.add(lineValues);
            }
        }

        processPair();

        System.out.println();
        System.out.println();

        System.out.println(output.toString());
    }

    private static final List<String> featureToIgnore = Stream.of("plip_pistack",
            "plip_cation",
            "plip_salt",
            "plip_metal",
            "plip_halogen",
            "plip_water",
            "plip_bb",
            "plip_total",
            "plip_l_hydrogen",
            "plip_l_hydrophobic",
            "plip_l_pistack",
            "plip_l_cation",
            "plip_l_salt",
            "plip_l_metal",
            "plip_halogen",
            "plip_water",
            "plip_l_bb",
            "plip_l_total",
            "plip_nl_hydrogen",
            "plip_nl_hydrophobic",
            "plip_nl_pistack",
            "plip_nl_cation",
            "plip_nl_salt",
            "plip_nl_metal",
            "plip_halogen",
            "plip_water",
            "plip_nl_bb",
            "plip_nl_total",
            "plip_betweenness",
            "plip_closeness",
            "plip_clusteringcoefficient",
            "plip_hydrogen_betweenness",
            "plip_hydrogen_closeness",
            "plip_hydrogen_clusteringcoefficient",
            "plip_hydrophobic_betweenness",
            "plip_hydrophobic_closeness",
            "plip_hydrophobic_clusteringcoefficient",
            "plip_distinct_neighborhoods")
            .collect(Collectors.toList());

    private static void processPair() {
        System.out.println("new entry, last pair was " + lastThermophilePdbId + " vs " + lastMesophilePdbId);
        String pairId = lastThermophilePdbId.toLowerCase() + "-" + lastMesophilePdbId.toLowerCase();

        StringJoiner pairOutput = new StringJoiner(System.lineSeparator());
        boolean[] sane = new boolean[] { true };

        IntStream.range(0, featureNames.size())
                .forEach(i -> {
                    String featureName = featureNames.get(i);
                    if(featureToIgnore.contains(featureName)) {
                        return;
                    }

                    double thermophileAverage = thermophileProtein.stream()
                            .mapToDouble(list -> list.get(i))
                            .average()
                            .getAsDouble();
                    double mesophileAverage = mesophileProtein.stream()
                            .mapToDouble(list -> list.get(i))
                            .average()
                            .getAsDouble();
                    double change = thermophileAverage / mesophileAverage * 100 - 100;

                    if(thermophileAverage == 0 || mesophileAverage == 0) {
                        sane[0] = false;
                    }

                    System.out.println(featureName + ": " +
                            StandardFormat.format(thermophileAverage) + " " +
                            StandardFormat.format(mesophileAverage) + " " +
                            StandardFormat.formatToInteger(change) + "%");

                    pairOutput.add(pairId + "," + featureName + "," + StandardFormat.formatToInteger(change));
        });

        if(sane[0]) {
            output.add(pairOutput.toString());
        }

        thermophileProtein.clear();
        mesophileProtein.clear();
    }
}
