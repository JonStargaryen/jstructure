package aars;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

/**
 * Path constants for the aars project.
 * Created by bittrich on 2/7/17.
 */
@Deprecated
class AARSConstants {
    static final String HOME_PATH = System.getProperty("user.home");
    private static final String PROJECT_BASE_PATH = HOME_PATH + "/git/aars_analysis/data/";
    static final String ARGTWEEZER_GEOMETRY_PATH = PROJECT_BASE_PATH + "geometry/argtweezer_geometry.tsv";
    static final String RENUMBERED_STRUCTURES_PATH = PROJECT_BASE_PATH + "msa/";
    static final String RENUMBERED_STRUCTURES_C2_PATH = RENUMBERED_STRUCTURES_PATH + "C2/renumbered_structures/";
    static final String MAIN_TABLE_CATALYTIC_PATH = PROJECT_BASE_PATH + "aars_main_table_catalytic.csv";
    static final String BINDING_SITE_PATH = PROJECT_BASE_PATH + "binding_sites/";

    static Stream<String> lines(Path path) {
        try {
            return Files.lines(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static Stream<Path> list(Path path) {
        try {
            return Files.list(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, byte[] bytes) {
        try {
            Files.createDirectories(path.getParent());
            Files.write(path, bytes);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
