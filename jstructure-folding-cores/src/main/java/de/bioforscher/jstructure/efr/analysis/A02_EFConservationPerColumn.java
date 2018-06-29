package de.bioforscher.jstructure.efr.analysis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Compute the EF conservation of each column in the MSA. Are some positions highly EFing or are the predictions
 * randomly distributed along the sequence?
 * best decision thresholds: 0.163 for the estimated probabilities.
 */
public class A02_EFConservationPerColumn {
    public static void main(String[] args) throws IOException {
        List<Double> class1Conservation = extractSequenceConservationFromFile("/home/bittrich/git/aars_data/T05_msa/C1_conservation.txt");
        List<Double> class2Conservation = extractSequenceConservationFromFile("/home/bittrich/git/aars_data/T05_msa/C2_conservation.txt");

        System.out.println(class1Conservation.size());
        System.out.println(class2Conservation.size());
    }

    private static List<Double> extractSequenceConservationFromFile(String filename) throws IOException {
        String line = Files.lines(Paths.get(filename))
                .filter(l -> l.startsWith("BAR_GRAPH\tConservation\tConservation of total alignment less than 25% gaps\t"))
                .findFirst()
                .get()
                .split("\t")[3];

        return Pattern.compile("\\|").splitAsStream(line)
                .map(entry -> entry.split(",")[0])
                .map(Double::valueOf)
                .collect(Collectors.toList());
    }
}
