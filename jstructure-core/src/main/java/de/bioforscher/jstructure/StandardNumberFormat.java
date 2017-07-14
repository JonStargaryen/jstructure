package de.bioforscher.jstructure;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/**
 * Project-wide way to format decimal numbers uniquely.
 * Created by bittrich on 7/13/17.
 */
public class StandardNumberFormat {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.US));

    public static DecimalFormat getInstance() {
        return DECIMAL_FORMAT;
    }

    public static String format(double value) {
        return DECIMAL_FORMAT.format(value);
    }

    public static String format(long value) {
        return DECIMAL_FORMAT.format(value);
    }

    public static String format(Object object) {
        return DECIMAL_FORMAT.format(object);
    }
}
