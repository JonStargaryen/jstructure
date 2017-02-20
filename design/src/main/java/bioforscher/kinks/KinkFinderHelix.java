package bioforscher.kinks;

/**
 * Representation of KinkFinder.
 * Created by bittrich on 2/14/17.
 */
public class KinkFinderHelix {
    private final String pdbCode, sequence;
    private final int helixStart, helixEnd, kinkPosition, kinkStart, kinkEnd;
    private final double kinkAngle, nRadius, nRmsd, cRadius, cRmsd, estimatedError;

    public KinkFinderHelix(String pdbCode, int helixStart, int helixEnd, int kinkPosition, int kinkStart, int kinkEnd, double kinkAngle, String sequence, double nRadius, double nRmsd, double cRadius, double cRmsd, double estimatedError) {
        this.pdbCode = pdbCode;
        this.helixStart = helixStart;
        this.helixEnd = helixEnd;
        this.kinkPosition = kinkPosition;
        this.kinkStart = kinkStart;
        this.kinkEnd = kinkEnd;
        this.kinkAngle = kinkAngle;
        this.sequence = sequence;
        this.nRadius = nRadius;
        this.nRmsd = nRmsd;
        this.cRadius = cRadius;
        this.cRmsd = cRmsd;
        this.estimatedError = estimatedError;
    }

    KinkFinderHelix(String[] lineToParse) {
        this(lineToParse[0],
                Integer.valueOf(lineToParse[1]),
                Integer.valueOf(lineToParse[2]),
                Integer.valueOf(lineToParse[3]),
                Integer.valueOf(lineToParse[4]),
                Integer.valueOf(lineToParse[5]),
                Double.valueOf(lineToParse[6]),
                lineToParse[7],
                Double.valueOf(lineToParse[8]),
                Double.valueOf(lineToParse[9]),
                Double.valueOf(lineToParse[10]),
                Double.valueOf(lineToParse[11]),
                Double.valueOf(lineToParse[12]));
    }

    public int getHelixStart() {
        return helixStart;
    }

    public int getHelixEnd() {
        return helixEnd;
    }

    public int getKinkPosition() {
        return kinkPosition;
    }

    public int getKinkStart() {
        return kinkStart;
    }

    public int getKinkEnd() {
        return kinkEnd;
    }

    public double getKinkAngle() {
        return kinkAngle;
    }

    public double getnRadius() {
        return nRadius;
    }

    public double getnRmsd() {
        return nRmsd;
    }

    public double getcRadius() {
        return cRadius;
    }

    public double getcRmsd() {
        return cRmsd;
    }

    public double getEstimatedError() {
        return estimatedError;
    }

    public String getSequence() {
        return sequence;
    }

    public String getPdbCode() {
        return pdbCode;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() +
                " sequence='" + sequence + '\'' +
                " helixStart=" + helixStart +
                " helixEnd=" + helixEnd +
                " kinkPosition=" + kinkPosition +
                " kinkAngle=" + kinkAngle;
    }

    public boolean isSignificantKink() {
        return kinkAngle > 20.0;
    }
}
