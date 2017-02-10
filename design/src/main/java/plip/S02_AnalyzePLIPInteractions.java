package plip;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Analyze retrieved PLIP data.
 * Created by bittrich on 2/10/17.
 */
public class S02_AnalyzePLIPInteractions {
    public static void main(String[] args) throws IOException {
        Files.walk(Paths.get(S01_MinePLIPInteractions.BASE_PATH + "interactions/"))
                .filter(path -> path.toFile().isFile())
                .forEach(System.out::println);
    }
}
