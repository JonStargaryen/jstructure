package bioforscher;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.function.Predicate;
import java.util.stream.Stream;

/**
 * Project constants.
 * Created by bittrich on 2/13/17.
 */
public class Constants {
    public static final String USER_HOME = System.getProperty("user.home");
    public static final String PROJECT_PATH = USER_HOME + "/git/phd_sb_repo/";
    public static final String DATA_PATH = PROJECT_PATH + "data/";
    public static final String STRUCTURE_PATH = DATA_PATH + "structures/";
    public static final String CHAIN_PATH = DATA_PATH + "chains/";
    public static final String PDBTM_ALPHA_NR_LIST =  DATA_PATH + "pdbtm_alpha_nr.list.txt";
    public static final String SCRIPTS_PATH = DATA_PATH + "scripts/";
    public static final String KINK_FINDER_SCRIPT_PATH = SCRIPTS_PATH + "kink_finder/Kink_Finder.py";
    public static final String KINK_FINDER_RESULT_PATH = DATA_PATH + "kink_finder/";
    public static final String OPM_PATH = DATA_PATH + "opm/";
    public static final String KINK_FINDER_TMP_PATH = KINK_FINDER_RESULT_PATH + "tmp/";
    public static final String STATISTICS_PATH = DATA_PATH + "statistics/";

    /**
     * The PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go).
     */
    public static final String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";

    /**
     * The OPM URL which can be used to fetch information by ID (format this using the id and you are good to go).
     */
    public static final String OPM_FETCH_URL = "http://opm.phar.umich.edu/protein.php?pdbid=%s";

    /**
     * Marks a line in a file as comment.
     */
    private static final String COMMENT_PREFIX = "#";

    /**
     * The suffix of PDB files.
     */
    public static final String PDB_SUFFIX = ".pdb";

    /**
     * True iff a given line starts with the comment prefix.
     */
    public static final Predicate<String> isCommentLine = line -> line.startsWith(COMMENT_PREFIX);

    public static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.US));

    private static void makeDirectoryIfAbsent(Path path) {
        try {
            if(!Files.exists(path)) {
                Files.createDirectories(path);
            }
        } catch (IOException e) {
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

    public static final String DELIMITER = "\t";

    public static Stream<Path> list(Path dir) {
        try {
            return Files.list(dir);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, byte[] bytes) {
        try {
            makeDirectoryIfAbsent(path.getParent());
            Files.write(path, bytes, StandardOpenOption.CREATE);
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
}
