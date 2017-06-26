package studies.gmlvq.fingerprint;

import studies.StudyConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Creates new version of itemset miner data set for GMLVQ-project.
 * Created by bittrich on 6/26/17.
 */
public class S01_ConvertFingerprintDataToArff {
    private static final Path basePath = StudyConstants.GMLVQ_MAIN
            .resolve("data")
            .resolve("fingerprint_miner")
            .resolve("results");

    public static void main(String[] args) throws IOException {
        // traverse over all directories in the base path - assumes identically structures subdirectories
        Files.list(basePath)
                .peek(path -> System.out.println("handling " + path))
                .flatMap(StudyConstants::list)
                .filter(Files::isDirectory)
                .peek(path -> System.out.println("handling " + path))
                .forEach(FingerPrintDataSetComposer::new);
    }
}
