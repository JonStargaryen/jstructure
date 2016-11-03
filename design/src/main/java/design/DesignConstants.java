package design;

import java.util.function.Predicate;

/**
 * Collection of paths and constants.
 * Created by S on 29.10.2016.
 */
public class DesignConstants {
    // directories
    /**
     * The root directory containing all data.
     */
    public static final String BASE_DIR = "D:/membrane/";

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

    /**
     * The directory containing all extracted sequences.
     */
    public static final String EXTRACTED_SEQUENCES_DIR = BASE_DIR + "extracted-sequences/";

    /**
     * All extracted sequences by topology.
     */
    public static final String EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR = EXTRACTED_SEQUENCES_DIR + "by-topology/";

    /**
     * The directory containing data as *.arff file.
     */
    public static final String ARFF_DIR = BASE_DIR + "arff/";

    /**
     * Where the sequences grouped by length live.
     */
    public static final String MOTIF_SEQUENCES_BY_LENGTH_DIR = ARFF_DIR + "motif-sequences-by-lengths/";

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

    // general predicates
    /**
     * True iff a given line starts with the comment prefix.
     */
    public static final Predicate<String> isCommentLine = line -> line.startsWith(COMMENT_PREFIX);
}
