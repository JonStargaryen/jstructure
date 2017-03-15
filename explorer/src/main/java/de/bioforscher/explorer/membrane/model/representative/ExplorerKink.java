package de.bioforscher.explorer.membrane.model.representative;

import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderHelix;

/**
 * The representation of kinks.
 * Created by bittrich on 3/1/17.
 */
public class ExplorerKink {
    private String chain, sequence;
    private int start, kinkPosition, end;
    private double angle, error;

    public ExplorerKink() {

    }

    ExplorerKink(KinkFinderHelix kinkFinderHelix) {
        this(kinkFinderHelix.getPdbCode().substring(4),
                kinkFinderHelix.getSequence(),
                kinkFinderHelix.getKinkPosition(),
                kinkFinderHelix.getKinkStart(),
                kinkFinderHelix.getKinkEnd(),
                kinkFinderHelix.getKinkAngle(),
                kinkFinderHelix.getEstimatedError());
    }

    ExplorerKink(String chain, String sequence, int start, int kinkPosition, int end, double angle, double error) {
        this.chain = chain;
        this.sequence = sequence;
        this.start = start;
        this.kinkPosition = kinkPosition;
        this.end = end;
        this.angle = angle;
        this.error = error;
    }

    public String getChain() {
        return chain;
    }

    public String getSequence() {
        return sequence;
    }

    public int getStart() {
        return start;
    }

    public int getKinkPosition() {
        return kinkPosition;
    }

    public int getEnd() {
        return end;
    }

    public double getAngle() {
        return angle;
    }

    public double getError() {
        return error;
    }
}
