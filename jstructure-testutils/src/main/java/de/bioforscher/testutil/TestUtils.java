package de.bioforscher.testutil;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
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

    public static Stream<String> getResourceAsStream(String filename) throws IOException {
        try(InputStreamReader inputStreamReader = new InputStreamReader(getResourceAsInputStream(filename))) {
            try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                return bufferedReader.lines();
            }
        }
    }

    public static List<String> getResourceAsLines(String filename) throws IOException {
        return getResourceAsStream(filename)
                .collect(Collectors.toList());
    }
}
