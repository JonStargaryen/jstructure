package de.bioforscher.jstructure.align.result.score;

abstract class AbstractAlignmentScore implements AlignmentScore {
    private final double score;

    AbstractAlignmentScore(double score) {
        this.score = score;
    }

    @Override
    public double getScore() {
        return score;
    }
}
