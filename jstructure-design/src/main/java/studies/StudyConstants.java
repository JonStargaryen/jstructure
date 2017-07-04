package studies;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

/**
 * Global constants and functions shared by all classes.
 * Created by bittrich on 5/17/17.
 */
public class StudyConstants {
    public static final Path HOME = Paths.get(System.getProperty("user.home"));
    public static final Path GIT = HOME.resolve("git");
    public static final Path GMLVQ_MAIN = GIT.resolve("gmlvq_main");
    public static final Path FINGERPRINT_MINER_PFAM_ARFF = StudyConstants.GMLVQ_MAIN.resolve("data").resolve("fingerprint_miner_pfam_arff");
    public static final Path PHD_REPO = GIT.resolve("phd_sb_repo");
    public static final Path AHAH_DATA_SET = PHD_REPO.resolve("data").resolve("AHAH").resolve("data").resolve("ahah_gold_standard.csv");

    public static Stream<Path> list(Path path) {
        try {
            return Files.list(path);
        } catch(IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<String> lines(Path path) {
        try {
            return Files.lines(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, String content) {
        try {
            Files.write(path, content.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void delete(Path path) {
        try {
            Files.delete(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void createDirectories(Path directory) {
        try {
            Files.createDirectories(directory);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<Path> walk(Path path) {
        try {
            return Files.walk(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
