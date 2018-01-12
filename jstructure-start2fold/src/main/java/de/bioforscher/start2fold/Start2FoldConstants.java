package de.bioforscher.start2fold;

import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;

public class Start2FoldConstants extends FileUtils {
    public static final String START2FOLD_URL = "http://start2fold.eu/";
    public static final String START2FOLD_ENTRIES_URL = "http://start2fold.eu/ids";
    public static final String START2FOLD_PROTEINS_URL = "http://start2fold.eu/proteins";

    public static final Path BASE_DIRECTORY = FileUtils.DATA_DIRECTORY.resolve("start2fold");
    public static final Path XML_DIRECTORY = BASE_DIRECTORY.resolve("xml");
    public static final Path FASTA_DIRECTORY = BASE_DIRECTORY.resolve("fasta");
    public static final Path PDB_DIRECTORY = BASE_DIRECTORY.resolve("pdb");
    public static final Path COUPLING_DIRECTORY = BASE_DIRECTORY.resolve("coupling");
    public static final Path PANCSA_LIST = BASE_DIRECTORY.resolve("pancsa.list");
    public static final Path PYMOL_DIRECTORY = BASE_DIRECTORY.resolve("pymol");
    public static final Path EQUANT_DIRECTORY = BASE_DIRECTORY.resolve("equant");
    public static final Path STATISTICS_DIRECTORY = BASE_DIRECTORY.resolve("statistics");
}
