package studies.gmlvq.csa;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Check if one CSA-entry has multiple EC-numbers.
 */
public class S02_ValidateDateset {
    public static void main(String[] args) throws IOException {
        // does one csa entry map to multiple ec numbers?
        Files.list(Paths.get("/home/bittrich/csa/2016-03-01_CSA_motifs/"))
                .filter(path -> !path.toFile().getName().equals("summary.csv"))
                .forEach(S02_ValidateDateset::handleFile);
    }

    private static void handleFile(Path path) {
        try {
            String id = path.toFile().getName().split("\\.")[0];
            if(Files.lines(path)
                    .filter(line -> line.startsWith("REMARK EC"))
                    .count() != 1) {
                System.err.println(id);
            }
        } catch (Exception e) {

        }
    }
}
