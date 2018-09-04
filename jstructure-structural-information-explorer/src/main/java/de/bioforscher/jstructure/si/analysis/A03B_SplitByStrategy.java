package de.bioforscher.jstructure.si.analysis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class A03B_SplitByStrategy {
    private static final Path BASE_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/statistics/");

    public static void main(String[] args) throws IOException {
        // create general bins
        String general = Files.lines(BASE_PATH.resolve("strategy.csv"))
                .filter(line -> !line.startsWith("id"))
                .filter(A03B_SplitByStrategy::isGeneralBin)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,strategy,rmsd" + System.lineSeparator(),
                        ""));
        Files.write(BASE_PATH.resolve("strategy-general.csv"),
                general.getBytes());

        // create best with false positives bins
        String best = Files.lines(BASE_PATH.resolve("strategy.csv"))
                .filter(line -> !line.startsWith("id"))
                .filter(A03B_SplitByStrategy::isBestFpBin)
                .map(A03B_SplitByStrategy::transformBestFpBin)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,fpperc,rmsd" + System.lineSeparator(),
                        ""));
        Files.write(BASE_PATH.resolve("strategy-bestfp.csv"),
                best.getBytes());

        // create random with false positive bins
        String random = Files.lines(BASE_PATH.resolve("strategy.csv"))
                .filter(line -> !line.startsWith("id"))
                .filter(A03B_SplitByStrategy::isRandomBin)
                .map(A03B_SplitByStrategy::transformRandomFpBin)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,fpperc,rmsd" + System.lineSeparator(),
                        ""));
        Files.write(BASE_PATH.resolve("strategy-randomfp.csv"),
                random.getBytes());

        // create worst with false positive bins
        String worst = Files.lines(BASE_PATH.resolve("strategy.csv"))
                .filter(line -> !line.startsWith("id"))
                .filter(A03B_SplitByStrategy::isWorstFpBin)
                .map(A03B_SplitByStrategy::transformWorstFpBin)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,fpperc,rmsd" + System.lineSeparator(),
                        ""));
        Files.write(BASE_PATH.resolve("strategy-worstfp.csv"),
                worst.getBytes());
    }

    private static String transformBestFpBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        String perc;
        if(strategy.equals("best")) {
            perc = "0";
        } else {
            perc = strategy.split("nonnative")[1];
        }
        return split[0] + "," + perc + "," + split[2];
    }

    private static String transformWorstFpBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        String perc;
        if(strategy.equals("worst")) {
            perc = "0";
        } else {
            perc = strategy.split("nonnative")[1];
        }
        return split[0] + "," + perc + "," + split[2];
    }

    private static String transformRandomFpBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        String perc;
        if(strategy.equals("random")) {
            perc = "0";
        } else {
            perc = strategy.split("nonnative")[1];
        }
        return split[0] + "," + perc + "," + split[2];
    }

    private static boolean isBestFpBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        return strategy.startsWith("best");
    }

    private static boolean isWorstFpBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        return strategy.startsWith("worst") || strategy.equals("best0_nonnative100");
    }

    private static boolean isRandomBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        return strategy.startsWith("random") || strategy.equals("best0_nonnative100");
    }

    private static boolean isGeneralBin(String line) {
        String[] split = line.split(",");
        String strategy = split[1];
        return strategy.equals("random") ||
                strategy.contains("short") ||
                strategy.contains("long") ||
                strategy.equals("best") ||
                strategy.equals("worst") ||
                strategy.contains("hydrogen") ||
                strategy.contains("hydrophobic");
    }
}
