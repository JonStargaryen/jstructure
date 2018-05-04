package de.bioforscher.jstructure.efr.collect.si;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import org.apache.commons.math3.stat.StatUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class A03_ReportOfRemovalAndAddition {
    private static List<Double> removal;
    private static List<Double> addition;

    public static void main(String[] args) throws IOException {
        removal = new ArrayList<>();
        addition = new ArrayList<>();

        Files.lines(Start2FoldConstants.BASE_DIRECTORY.resolve("pancsa-si.list"))
                .forEach(A03_ReportOfRemovalAndAddition::handleLine);

        double[] r = removal.stream()
                .mapToDouble(Double::valueOf)
                .toArray();
        double[] a = addition.stream()
                .mapToDouble(Double::valueOf)
                .toArray();

        System.out.println("removal");
        System.out.println(StatUtils.mean(r));
        System.out.println(Math.sqrt(StatUtils.variance(r)));

        System.out.println("addition");
        System.out.println(StatUtils.mean(a));
        System.out.println(Math.sqrt(StatUtils.variance(a)));
    }

    private static void handleLine(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];

            Files.lines(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("raw").resolve(entryId.toUpperCase() + ".out"))
                    .forEach(l -> {
                        boolean wasRemoval = l.contains("true");
                        double rmsd = Double.valueOf(l.split("\t")[8]);
                        if(wasRemoval) {
                            removal.add(rmsd);
                        } else {
                            addition.add(rmsd);
                        }
                    });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
