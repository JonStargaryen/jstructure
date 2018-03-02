package de.bioforscher.jstructure;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.Locale;

/**
 * Project-wide way to format decimal numbers uniquely.
 * Created by bittrich on 7/13/17.
 */
public class StandardFormat {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.US));
    private static final DecimalFormat INTEGER_FORMAT = new DecimalFormat("#", new DecimalFormatSymbols(Locale.US));

//    private static DateTimeFormatter PDB_DATE_FORMAT = DateTimeFormatter.ofPattern("dd-MMM-yy");
private static DateTimeFormatter PDB_DATE_FORMAT = new DateTimeFormatterBuilder()
        .parseCaseInsensitive()
        .appendPattern("dd-MMM-yy")
        .toFormatter(Locale.US);

    public static DecimalFormat getDecimalFormatInstance() {
        return DECIMAL_FORMAT;
    }

    public static DateTimeFormatter getPdbDateFormatInstance() {
        return PDB_DATE_FORMAT;
    }

    public static String format(double value) {
        return DECIMAL_FORMAT.format(value);
    }

    public static String formatToInteger(double value) {
        return INTEGER_FORMAT.format(value);
    }

    public static String format(long value) {
        return DECIMAL_FORMAT.format(value);
    }

    public static String format(Object object) {
        return DECIMAL_FORMAT.format(object);
    }

    public static String format(LocalDate localDate) {
        return localDate.format(PDB_DATE_FORMAT).toUpperCase();
    }
}
