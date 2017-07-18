package de.bioforscher.jstructure.mutation.impl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Wraps the rate4site executable (which is wrapped in a Docker container).
 * Created by bittrich on 7/17/17.
 */
public class Rate4SiteWrapper {
    List<Double> executeCommand(String alignmentString) throws IOException, InterruptedException {
        // write to tmp file
        Path alignmentFile = Files.createTempFile("r4s", "aln");
        Files.write(alignmentFile, alignmentString.getBytes());

        // execute command
        ProcessBuilder processBuilder = new ProcessBuilder("docker",
                "run",
                "--rm",
                "-v",
                "/tmp/:/mnt",
                "bigm/rate4site",
                "-s",
                "/mnt/" + alignmentFile.toFile().getName());
        Process process = processBuilder.start();

        List<Double> output = new BufferedReader(new InputStreamReader(process.getInputStream()))
                .lines()
                // drop additional output lines
                .filter(line -> line.contains("rate of pos"))
                .map(line -> line.split("=")[1])
                .map(String::trim)
                .map(Double::valueOf)
                .collect(Collectors.toList());
        process.waitFor();

        // delete tmp file
        Files.delete(alignmentFile);

        return output;
    }
}
