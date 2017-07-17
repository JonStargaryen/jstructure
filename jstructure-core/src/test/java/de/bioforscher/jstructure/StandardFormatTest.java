package de.bioforscher.jstructure;

import org.junit.Assert;
import org.junit.Test;

import java.time.LocalDate;

/**
 * Tests for the project-wide formatting functions.
 * Created by S on 14.07.2017.
 */
public class StandardFormatTest {
    @Test
    public void shouldFormatDateAccordingToPdbHeaderFormat() {
        LocalDate date = LocalDate.now();
        System.out.println(StandardFormat.format(date));
    }

    @Test
    public void shouldParseToPreviousCentury() {
        LocalDate date = LocalDate.parse("01-JAN-93", StandardFormat.getPdbDateFormatInstance());
        if(date.isAfter(LocalDate.now())) {
            date = date.minusYears(100);
        }
        System.out.println(date);
        Assert.assertTrue("LocalDate of e.g. 01-JAN-93 is parsed to 2093-01-01",
                date.isBefore(LocalDate.now()));
    }
}