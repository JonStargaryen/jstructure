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
        Assert.fail("LocalDate of e.g. 93-JAN-01 is parsed to 2093-01-01");
    }
}