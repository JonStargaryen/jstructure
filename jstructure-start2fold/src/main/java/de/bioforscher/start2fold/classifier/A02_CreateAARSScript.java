package de.bioforscher.start2fold.classifier;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

public class A02_CreateAARSScript {
    public static void main(String[] args) throws IOException {
        String output = Files.lines(Paths.get("/home/bittrich/git/aars_data/dataset.csv"))
                .filter(line -> !line.startsWith("identifier"))
                .map(line -> line.split(",")[0])
                .map(id -> id.split("_")[0])
                .distinct()
                .map(pdbId -> "java -jar efr-classifier.jar " + pdbId + " " + pdbId + ".out")
                .collect(Collectors.joining(System.lineSeparator()));

        Files.write(Paths.get("/home/bittrich/classifier/run.sh"),
                output.getBytes());
    }
}
