package de.bioforscher.start2fold.model;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

public class EQuantScore extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
    private final double evaluation;

    public EQuantScore(double evaluation) {
        super(null);
        this.evaluation = evaluation;
    }

    public double getEvaluation() {
        return evaluation;
    }

    @Override
    public Double getValue() {
        return getEvaluation();
    }
}
