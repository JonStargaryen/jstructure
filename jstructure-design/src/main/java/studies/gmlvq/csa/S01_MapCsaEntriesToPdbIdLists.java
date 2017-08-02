package studies.gmlvq.csa;

import de.bioforscher.jstructure.feature.mapping.SiftsMappingAnnotator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

public class S01_MapCsaEntriesToPdbIdLists {
    public static void main(String[] args) throws IOException {
        Files.list(Paths.get("/home/bittrich/csa/2016-03-01_CSA_motifs/"))
                .filter(path -> !path.toFile().getName().equals("summary.csv"))
                .forEach(S01_MapCsaEntriesToPdbIdLists::handleFile);
    }

    private static void handleFile(Path path) {
        try {
            String id = path.toFile().getName().split("\\.")[0];
            String ecNumber = Files.lines(path)
                    .filter(line -> line.startsWith("REMARK EC"))
                    .findFirst()
                    .get()
                    .split("EC")[1];
            List<String> pdbIds = SiftsMappingAnnotator.getLinesForEc(ecNumber)
                    .map(split -> split[0].toLowerCase() + "_" + split[1])
                    .collect(Collectors.toList());
            System.out.println(id + " " + ecNumber + ":\n" + pdbIds);
            Files.write(Paths.get("/home/bittrich/csa/pdb_mapping/").resolve(id + ".list"),
                    (ecNumber + System.lineSeparator() + pdbIds).getBytes());
        } catch (Exception e) {

        }
    }
}
