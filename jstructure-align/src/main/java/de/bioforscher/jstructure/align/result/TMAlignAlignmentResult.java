package de.bioforscher.jstructure.align.result;

import de.bioforscher.jstructure.align.result.score.RootMeanSquareDeviation;
import de.bioforscher.jstructure.align.result.score.TemplateModelingScore;

public class TMAlignAlignmentResult implements AlignmentResult {
    private final int length1;
    private final int length2;
    private final int alignedLength;
    private final RootMeanSquareDeviation rootMeanSquareDeviation;
    private final double seqId;
    private final TemplateModelingScore templateModelingScore1;
    private final TemplateModelingScore templateModelingScore2;

    public TMAlignAlignmentResult(int length1,
                                  int length2,
                                  int alignedLength,
                                  RootMeanSquareDeviation rootMeanSquareDeviation,
                                  double seqId,
                                  TemplateModelingScore templateModelingScore1,
                                  TemplateModelingScore templateModelingScore2) {
        this.length1 = length1;
        this.length2 = length2;
        this.alignedLength = alignedLength;
        this.rootMeanSquareDeviation = rootMeanSquareDeviation;
        this.seqId = seqId;
        this.templateModelingScore1 = templateModelingScore1;
        this.templateModelingScore2 = templateModelingScore2;
    }

    public int getLength1() {
        return length1;
    }

    public int getLength2() {
        return length2;
    }

    public int getAlignedLength() {
        return alignedLength;
    }

    public double getSeqId() {
        return seqId;
    }

    public RootMeanSquareDeviation getRootMeanSquareDeviation() {
        return rootMeanSquareDeviation;
    }

    public TemplateModelingScore getAverageTemplateModelingScore() {
        return new TemplateModelingScore(0.5 * templateModelingScore1.getScore() + 0.5 * templateModelingScore2.getScore());
    }

    public TemplateModelingScore getTemplateModelingScore1() {
        return templateModelingScore1;
    }

    public TemplateModelingScore getTemplateModelingScore2() {
        return templateModelingScore2;
    }

    @Override
    public String toString() {
        return "TMAlignAlignmentResult{" +
                "length1=" + length1 +
                ", length2=" + length2 +
                ", alignedLength=" + alignedLength +
                ", rootMeanSquareDeviation=" + rootMeanSquareDeviation +
                ", seqId=" + seqId +
                ", templateModelingScore1=" + templateModelingScore1 +
                ", templateModelingScore2=" + templateModelingScore2 +
                '}';
    }
}
