package studies.gmlvq.fingerprint;

import studies.StudyConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Annotate data set as pipeline.
 * Created by bittrich on 6/29/17.
 */
public class S02_FingerprintDataAnnotationPipeline {
    private static final Path basePath = StudyConstants.HOME
            .resolve("fingerprint_miner_out");

    public static void main(String[] args) throws IOException {
        // traverse over all directories in the base path - assumes identically structured subdirectories
        Files.list(basePath)
                .peek(path -> System.out.println("handling " + path))
                .flatMap(StudyConstants::list)
                .filter(Files::isDirectory)
                .filter(path -> !path.toFile().getName().equals("itemset-miner"))
                .peek(path -> System.out.println("handling " + path))
                .forEach(FingerPrintDataSetComposer::new);
    }
}
