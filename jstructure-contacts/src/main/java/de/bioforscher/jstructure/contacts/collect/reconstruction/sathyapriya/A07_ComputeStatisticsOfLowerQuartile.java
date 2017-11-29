package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.LowestQuartileSelector;

public class A07_ComputeStatisticsOfLowerQuartile {
    public static void main(String[] args) {
        new LowestQuartileSelector(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("reconstruction-performance.csv"));
    }
}
