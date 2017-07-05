package studies.gmlvq.fingerprint;

import studies.StudyConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static studies.StudyConstants.FINGERPRINT_MINER_PFAM_ARFF;

/**
 * Annotate data set as pipeline.
 * Created by bittrich on 6/29/17.
 */
public class S01_CreateArffsForFingerprintData {
    public static final Path INPUT_PATH = StudyConstants.HOME.resolve("fingerprint-miner");

    public static void main(String[] args) throws IOException {
        // traverse over all directories in the base path - assumes identically structured subdirectories
        Files.list(INPUT_PATH)
                .filter(Files::isDirectory)
                .filter(path -> !path.toFile().getName().equals("results"))
                .forEach(path -> new FingerPrintDataSetComposer(path, FINGERPRINT_MINER_PFAM_ARFF));
    }
}
