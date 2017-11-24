package de.bioforscher.equant;

import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;

public class EquantConstants extends FileUtils {
    public static final Path EQUANT_DATA_DIRECTORY = DATA_DIRECTORY.resolve("equant");
    public static final Path CASP_DATA_DIRECTORY = EQUANT_DATA_DIRECTORY.resolve("casp");
    public static final Path CASP9_DIRECTORY = CASP_DATA_DIRECTORY.resolve("CASP9");
}
