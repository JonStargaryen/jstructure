package de.bioforscher.jstructure.align.result.score;

import de.bioforscher.jstructure.StandardFormat;

public class RootMeanSquareDeviation extends AbstractAlignmentScore {
    public RootMeanSquareDeviation(double rootMeanSquareDeviation) {
        super(rootMeanSquareDeviation);
    }

    @Override
    public String toString() {
        return StandardFormat.format(getScore());
    }
}
