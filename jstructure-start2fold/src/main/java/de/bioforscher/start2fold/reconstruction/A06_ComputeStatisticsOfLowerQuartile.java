package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;

public class A06_ComputeStatisticsOfLowerQuartile {
    public static void main(String[] args) {
        new LowestQuartileSelector(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("reconstruction-performance.csv"));
    }
}
