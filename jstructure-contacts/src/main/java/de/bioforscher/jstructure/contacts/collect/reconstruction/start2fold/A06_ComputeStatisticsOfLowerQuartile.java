package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.LowestQuartileSelector;

public class A06_ComputeStatisticsOfLowerQuartile {
    public static void main(String[] args) {
        new LowestQuartileSelector(ContactsConstants.START2FOLD_DIRECTORY.resolve("reconstruction-performance.csv"));
    }
}
