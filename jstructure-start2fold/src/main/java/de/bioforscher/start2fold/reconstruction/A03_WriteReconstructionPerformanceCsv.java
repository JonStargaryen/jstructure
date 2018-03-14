package de.bioforscher.start2fold.reconstruction;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A03_WriteReconstructionPerformanceCsv {
    private static final Path RECONSTRUCTION_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/");

    public static void main(String[] args) throws IOException {
        String output = Files.list(RECONSTRUCTION_PATH)
                .filter(Files::isRegularFile)
                .filter(path -> path.toFile().getName().endsWith(".out"))
                .map(A03_WriteReconstructionPerformanceCsv::handleFile)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id,strategy,percentage,mapId,modelId,rmsd" + System.lineSeparator(),
                        ""));

        Files.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/statistics/reconstruction.csv"),
                output.getBytes());
    }

    private static String handleFile(Path path) {
        try {
            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
            List<String> lines = Files.readAllLines(path);

            String idLine = "";
            for(String line : lines) {
                if(!line.startsWith("Aligned")) {
                    idLine = line;
                    continue;
                } else {
                    String rmsd = line.split("=")[2].split(",")[0].trim();
                    stringJoiner.add(idLine + "," + rmsd);
                }
            }

            return stringJoiner.toString();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
