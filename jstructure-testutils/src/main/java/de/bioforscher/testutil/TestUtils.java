package de.bioforscher.testutil;

import java.io.*;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.Path;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Access to commonly used test functions and constants.
 * Created by S on 29.09.2016.
 */
public class TestUtils {
    public static final double[] ZERO_VECTOR = new double[] { 0, 0, 0 };
    public static final double TOLERANT_ERROR_MARGIN = 0.001;

    public enum SupportedProtein {
        PDB_1ACJ,
        PDB_1AR1,
        PDB_1AEY,
        PDB_1BRR,
        PDB_1CYO,
        PDB_2LZM,
        PDB_5OAZ
    }

    /**
     * Fast access to protein instance to be used by tests.
     * Usage:
     * <code>Protein protein = ProteinParser.source(getProteinInputStream(PDB_1BRR))
     *                                      .minimalParsing(true)
     *                                      .parse()</code>
     * @param supportedProtein the protein's enum entry
     * @return an input stream of that pdb file's content
     */
    public static InputStream getProteinInputStream(SupportedProtein supportedProtein) {
        return getResourceAsInputStream("pdb/" + supportedProtein.name().split("_")[1] + ".pdb");
    }

    public static InputStream getResourceAsInputStream(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        InputStream is = ccl.getResourceAsStream(filename);
        return Objects.requireNonNull(is, "could not acquire inputstream of resource '" + filename + "'");
    }

    public static Stream<String> getResourceAsStream(String filename) {
        return getResourceAsLines(filename).stream();
    }

    public static List<String> getResourceAsLines(String filename) {
        try {
            try (InputStreamReader inputStreamReader = new InputStreamReader(getResourceAsInputStream(filename))) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    return bufferedReader.lines()
                            .collect(Collectors.toList());
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Path downloadPdbId(Path targetDirectory, String pdbId) {
        try {
            Path targetPath = targetDirectory.resolve(pdbId + ".pdb");
            URL url = new URL("https://files.rcsb.org/download/" + pdbId.toUpperCase() + ".pdb");
            ReadableByteChannel rbc = Channels.newChannel(url.openStream());
            FileOutputStream fos = new FileOutputStream(targetPath.toFile());
            fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
            return targetPath;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
