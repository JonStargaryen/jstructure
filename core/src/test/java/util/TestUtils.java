package util;

import java.net.URL;
import java.util.Objects;

/**
 * Created by S on 29.09.2016.
 */
public class TestUtils {
    public static String getResourceAsFilepath(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        URL resource = ccl.getResource(filename);
        Objects.requireNonNull(resource);
        return resource.getPath();
    }
}
