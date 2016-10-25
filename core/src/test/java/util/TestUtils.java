package util;

import java.net.URL;
import java.util.Objects;

/**
 * Access to commonly used test functions and constants.
 * Created by S on 29.09.2016.
 */
public class TestUtils {
    public static final double[] ZERO_VECTOR = new double[] { 0, 0, 0 };
    public static final double TOLERANT_ERROR_MARGIN = 0.001;

    public static String getResourceAsFilepath(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        URL resource = ccl.getResource(filename);
        Objects.requireNonNull(resource);
        return resource.getPath();
    }
}
