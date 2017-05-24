package studies;

import java.net.URL;
import java.util.Objects;

/**
 * Global constants and functions shared by all classes.
 * Created by bittrich on 5/17/17.
 */
public class StudyConstants {
    public static final String HOME = System.getProperty("user.home");
    public static final String GIT = HOME + "/git/";
    public static final String GMLVQ_MAIN = GIT + "gmlvq_main/";
    public static final String PHD = GIT + "phd_sb_repo/";

    public static String getResourceAsFilepath(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        URL resource = ccl.getResource(filename);
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Objects.requireNonNull(resource).getPath().replaceFirst("^/(.:/)", "$1");
    }
}
