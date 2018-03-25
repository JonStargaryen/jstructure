package de.bioforscher.jstructure.align.result;

import de.bioforscher.jstructure.align.result.score.RootMeanSquareDeviation;
import de.bioforscher.jstructure.align.result.score.TemplateModelingScore;

public class TMAlignAlignmentResult implements AlignmentResult {
    private final RootMeanSquareDeviation rootMeanSquareDeviation;
    private final TemplateModelingScore templateModelingScore;

    public TMAlignAlignmentResult(RootMeanSquareDeviation rootMeanSquareDeviation,
                                  TemplateModelingScore templateModelingScore) {
        this.rootMeanSquareDeviation = rootMeanSquareDeviation;
        this.templateModelingScore = templateModelingScore;
    }

    public RootMeanSquareDeviation getRootMeanSquareDeviation() {
        return rootMeanSquareDeviation;
    }

    public TemplateModelingScore getTemplateModelingScore() {
        return templateModelingScore;
    }
}
