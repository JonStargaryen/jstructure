package de.bioforscher.jstructure.si.analysis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class A03A_CleanupOfWronglyComputedBins {
    private static final Path BASE_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/statistics/");

    public static void main(String[] args) throws IOException {
        // create general bins
        String general = Files.lines(BASE_PATH.resolve("strategy.csv"))
                .filter(line -> !line.startsWith("id"))
                .filter(line -> !isErroneousBin(line))
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,strategy,rmsd" + System.lineSeparator(),
                        ""));
        Files.write(BASE_PATH.resolve("strategy.csv"),
                general.getBytes());
    }

    private static boolean isErroneousBin(String line) {
        String strategy = line.split(",")[1];
        boolean isErroneousLine = strategy.equals("short") ||
                strategy.equals("long") ||
                strategy.equals("hydrogen") ||
                strategy.equals("hydrophobic");
        if(isErroneousLine) {
            System.out.println("dropping line: " + line);
        }
        return isErroneousLine;
    }
}
