package de.bioforscher.testutil;

import java.io.*;
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

    public static InputStream getResourceAsInputStream(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        InputStream is = ccl.getResourceAsStream(filename);
        return Objects.requireNonNull(is);
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
}
