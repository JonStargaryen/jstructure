package util;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

/**
 * Set log level globally.
 * Created by S on 08.11.2016.
 */
public class LoggingUtils {
    public LoggingUtils() {
        Logger.getRootLogger().setLevel(Level.ALL);
    }
}
