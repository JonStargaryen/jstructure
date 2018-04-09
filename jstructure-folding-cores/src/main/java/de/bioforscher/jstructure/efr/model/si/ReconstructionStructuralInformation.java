package de.bioforscher.jstructure.efr.model.si;

import de.bioforscher.jstructure.efr.model.ContactDistanceBin;

public class ReconstructionStructuralInformation {
    private final int residueIdentifier1;
    private final String aa1;
    private final int residueIdentifier2;
    private final String aa2;
    private final ContactDistanceBin contactDistanceBin;
    private final boolean contactWasRemoved;
    private final double baselineRmsd;
    private final double baselineTmScore;
    private final double baselineQ;
    private final double reconstructionRmsd;
    private final double reconstructionTmScore;
    private final double reconstructionQ;
    private final double rmsdIncrease;
    private final double tmScoreIncrease;
    private final double qIncrease;

    public ReconstructionStructuralInformation(int residueIdentifier1,
                                               String aa1,
                                               int residueIdentifier2,
                                               String aa2,
                                               ContactDistanceBin contactDistanceBin,
                                               boolean contactWasRemoved,
                                               double baselineRmsd,
                                               double baselineTmScore,
                                               double baselineQ,
                                               double reconstructionRmsd,
                                               double reconstructionTmScore,
                                               double reconstructionQ,
                                               double rmsdIncrease,
                                               double tmScoreIncrease,
                                               double qIncrease) {
        this.residueIdentifier1 = residueIdentifier1;
        this.aa1 = aa1;
        this.residueIdentifier2 = residueIdentifier2;
        this.aa2 = aa2;
        this.contactDistanceBin = contactDistanceBin;
        this.contactWasRemoved = contactWasRemoved;
        this.baselineRmsd = baselineRmsd;
        this.baselineTmScore = baselineTmScore;
        this.baselineQ = baselineQ;
        this.reconstructionRmsd = reconstructionRmsd;
        this.reconstructionTmScore = reconstructionTmScore;
        this.reconstructionQ = reconstructionQ;
        this.rmsdIncrease = rmsdIncrease;
        this.tmScoreIncrease = tmScoreIncrease;
        this.qIncrease = qIncrease;
    }

    public int getResidueIdentifier1() {
        return residueIdentifier1;
    }

    public String getAa1() {
        return aa1;
    }

    public int getResidueIdentifier2() {
        return residueIdentifier2;
    }

    public String getAa2() {
        return aa2;
    }

    public ContactDistanceBin getContactDistanceBin() {
        return contactDistanceBin;
    }

    public boolean isContactWasRemoved() {
        return contactWasRemoved;
    }

    public double getBaselineRmsd() {
        return baselineRmsd;
    }

    public double getBaselineTmScore() {
        return baselineTmScore;
    }

    public double getBaselineQ() {
        return baselineQ;
    }

    public double getReconstructionRmsd() {
        return reconstructionRmsd;
    }

    public double getReconstructionTmScore() {
        return reconstructionTmScore;
    }

    public double getReconstructionQ() {
        return reconstructionQ;
    }

    public double getRmsdIncrease() {
        return rmsdIncrease;
    }

    public double getTmScoreIncrease() {
        return tmScoreIncrease;
    }

    public double getqIncrease() {
        return qIncrease;
    }
}
