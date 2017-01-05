package design;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.function.Predicate;
import java.util.stream.Stream;

/**
 * Collection of paths and constants.
 * Created by S on 29.10.2016.
 */
public class DesignConstants {
    // directories
    /**
     * The root directory containing all data.
     */
    public static final String BASE_DIR = "/home/bittrich/membrane/";
//    public static final String BASE_DIR = "D:/membrane/";

    /**
     * The list of considered pdb ids of non-redundant alpha-helical membrane proteins.
     */
    public static final String NR_ALPHA_IDS = BASE_DIR + "pdbtm_alpha_nr.list";

    /**
     * The directory containing all pdb files.
     */
    public static final String PDB_DIR = BASE_DIR + "pdb/";

    /**
     * The directory of raw OPM data.
     */
    public static final String OPM_RAW_DIR = BASE_DIR + "opm-raw/";

    /**
     * Motif fragment directory.
     */
    public static final String MOTIF_FRAGMENT_DIR = BASE_DIR + "motif-fragments/";

    /**
     * Motif fragments grouped by opm topology.
     */
    public static final String MOTIF_FRAGMENT_BY_TOPOLOGY_DIR = MOTIF_FRAGMENT_DIR + "by-topology/";

    @Deprecated
    public static final String ALIGNED_MOTIF_FRAGMENT_DIR = BASE_DIR + "aligned-motif-fragments/";

    @Deprecated
    public static final String ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR = ALIGNED_MOTIF_FRAGMENT_DIR + "by-topology/";

    @Deprecated
    public static final String ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_SAMPLED_DIR = ALIGNED_MOTIF_FRAGMENT_DIR + "by-topology-sampled/";

    @Deprecated
    public static final String ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR = ALIGNED_MOTIF_FRAGMENT_DIR + "consensus/";

    /**
     * The directory containing all extracted sequences.
     */
    public static final String EXTRACTED_SEQUENCES_DIR = BASE_DIR + "extracted-sequences/";

    /**
     * All extracted sequences by topology.
     */
    public static final String EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR = EXTRACTED_SEQUENCES_DIR + "by-topology/";

    /**
     * All extracted sequences clustered by their respective structural cluster.
     */
    public static final String EXTRACTED_SEQUENCES_CONSENSUS_DIR = EXTRACTED_SEQUENCES_DIR + "consensus/";

    /**
     * The directory containing data as *.arff file.
     */
    public static final String ARFF_DIR = BASE_DIR + "arff/";

    /**
     * The directory containing statistics output.
     */
    public static final String STATISTICS_DIR = BASE_DIR + "statistics/";

    /**
     * Where the sequences grouped by length live.
     */
    public static final String MOTIF_SEQUENCES_BY_LENGTH_DIR = ARFF_DIR + "motif-sequences-by-lengths/";

    /**
     * The PyMOL runtime location.
     */
    public static final String PYMOL_DIR = "C:/Python27/PyMOL/";

    public static final String FRAGMENT_CLUSTERS = BASE_DIR + "fragment-clusters/";

    public static final String NAIVE_FRAGMENT_CLUSTERS = BASE_DIR + "naive-clusters/"   ;

    public static final String FRAGMENT_CONSENSUS_CLUSTERS = FRAGMENT_CLUSTERS + "consensus/";

    /**
     * The PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go).
     */
    public static final String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";

    /**
     * The OPM URL which can be used to fetch information by ID (format this using the id and you are good to go).
     */
    public static final String OPM_FETCH_URL = "http://opm.phar.umich.edu/protein.php?pdbid=%s";

    // file-related constants
    /**
     * Marks a line in a file as comment.
     */
    public static final String COMMENT_PREFIX = "#";

    /**
     * The suffix of PDB files.
     */
    public static final String PDB_SUFFIX = ".pdb";

    /**
     * The suffix of raw OPM files.
     */
    public static final String OPM_SUFFIX = ".opm";

    /**
     * The suffix of sequence files.
     */
    public static final String SEQUENCE_SUFFIX = ".seq";

    /**
     * The suffix of arff files.
     */
    public static final String ARFF_SUFFIX = ".arff";

    /**
     * The default delimiter.
     */
    public static final String DELIMITER = "\t";

    /**
     * The known topologies.
     */
    public static final List<String> TOPOLOGIES = Arrays.asList("tm", "ntm", "trans");

    /**
     * The threshold below a cluster is considered to be rare and will be merged with other seldom clusters to the
     * 'rare' cluster with id {@link DesignConstants#RARE_CLUSTER_NAME}.
     */
    public static final double RARE_CLUSTER_THRESHOLD = 0.05;

    /**
     * The threshold below 2 fragments are assumed to be of the same cluster.
     */
    public static final double IDENTICAL_CLUSTER_RMSD_THRESHOLD = 0.5;

    /**
     * The identifier of the merged cluster of seldom occurrences.
     */
    public static final String RARE_CLUSTER_NAME = "0";

    public static final String CLUSTER_SUMMARY = "clusters.summary";

    public static final String CLUSTER_CONSENSUS = "consensus" + PDB_SUFFIX;

    public static final String TSV_DIR = BASE_DIR + "tsv/";

    public static Stream<Path> list(Path dir) {
        try {
            return Files.list(dir);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, byte[] bytes) {
        try {
            Files.write(path, bytes, StandardOpenOption.CREATE);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    // general predicates
    /**
     * True iff a given line starts with the comment prefix.
     */
    public static final Predicate<String> isCommentLine = line -> line.startsWith(COMMENT_PREFIX);

    public static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.US));

    public static void makeDirectoryIfAbsent(Path path) {
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
}
