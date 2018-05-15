package de.bioforscher.jstructure.si.analysis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class A03A_FilterReconstructionStrategyCsv {
    public static void main(String[] args) throws IOException {
        String output = Files.lines(Paths.get("/home/bittrich/git/phd_sb_repo/data/si/statistics/reconstruction-strategy.csv"))
                .filter(line -> line.startsWith("id") || line.startsWith("STF"))
                .collect(Collectors.joining(System.lineSeparator()));

        Files.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/si/statistics/reconstruction-strategy-filtered.csv"),
                output.getBytes());
    }
}
