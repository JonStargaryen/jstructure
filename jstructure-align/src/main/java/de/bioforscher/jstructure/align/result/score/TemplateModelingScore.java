package de.bioforscher.jstructure.align.result.score;

import de.bioforscher.jstructure.StandardFormat;

public class TemplateModelingScore extends AbstractAlignmentScore {
    public TemplateModelingScore(double templateModelingScore) {
        super(templateModelingScore);
    }

    @Override
    public String toString() {
        return StandardFormat.format(getScore());
    }
}
