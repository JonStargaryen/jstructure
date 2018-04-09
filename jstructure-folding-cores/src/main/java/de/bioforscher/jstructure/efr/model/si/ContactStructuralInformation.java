package de.bioforscher.jstructure.efr.model.si;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.model.ContactDistanceBin;

import java.util.List;

public class ContactStructuralInformation {
    private final int residueIdentifier1;
    private final String aa1;
    private final int residueIdentifier2;
    private final String aa2;
    private final ContactDistanceBin contactDistanceBin;
    private final double baselineRmsd;
    private final double baselineTmScore;
    private final double baselineQ;
    private final double averageRmsdIncrease;
    private final double averageTmScoreIncrease;
    private final double averageQIncrease;
    private final double maximumRmsdIncrease;
    private final double maximumTmScoreIncrease;
    private final double maximumQIncrease;
    private final boolean isEarlyFoldingResidue;
    private final boolean isEarlyFoldingContact;
    private final double averageRmsdIncreaseZScore;
    private final double maximumRmsdIncreaseZScore;
    private final double fractionOfTopScoringContacts;
    private double plmScore;

    public ContactStructuralInformation(int residueIdentifier1,
                                        String aa1,
                                        int residueIdentifier2,
                                        String aa2,
                                        ContactDistanceBin contactDistanceBin,
                                        double baselineRmsd,
                                        double baselineTmScore,
                                        double baselineQ,
                                        double averageRmsdIncrease,
                                        double averageTmScoreIncrease,
                                        double averageQIncrease,
                                        double maximumRmsdIncrease,
                                        double maximumTmScoreIncrease,
                                        double maximumQIncrease,
                                        boolean isEarlyFoldingResidue,
                                        boolean isEarlyFoldingContact,
                                        double averageRmsd,
                                        double standardDeviationRmsd,
                                        double averageMaximumRmsd,
                                        double standardDeviationMaximumRmsd,
                                        List<ReconstructionStructuralInformation> allReconstructions,
                                        List<ReconstructionStructuralInformation> topScoringReconstructions) {
        this.residueIdentifier1 = residueIdentifier1;
        this.aa1 = aa1;
        this.residueIdentifier2 = residueIdentifier2;
        this.aa2 = aa2;
        this.contactDistanceBin = contactDistanceBin;
        this.baselineRmsd = baselineRmsd;
        this.baselineTmScore = baselineTmScore;
        this.baselineQ = baselineQ;
        this.averageRmsdIncrease = averageRmsdIncrease;
        this.averageTmScoreIncrease = averageTmScoreIncrease;
        this.averageQIncrease = averageQIncrease;
        this.maximumRmsdIncrease = maximumRmsdIncrease;
        this.maximumTmScoreIncrease = maximumTmScoreIncrease;
        this.maximumQIncrease = maximumQIncrease;
        this.isEarlyFoldingResidue = isEarlyFoldingResidue;
        this.isEarlyFoldingContact = isEarlyFoldingContact;
        this.averageRmsdIncreaseZScore = (averageRmsdIncrease - averageRmsd) / standardDeviationRmsd;
        this.maximumRmsdIncreaseZScore = (maximumRmsdIncrease - averageMaximumRmsd) / standardDeviationMaximumRmsd;
        double frac = topScoringReconstructions.stream()
                .filter(recon -> residueIdentifier1 == recon.getResidueIdentifier1() && residueIdentifier2 == recon.getResidueIdentifier2())
                .count() / (double) allReconstructions.stream()
                .filter(recon -> residueIdentifier1 == recon.getResidueIdentifier1() && residueIdentifier2 == recon.getResidueIdentifier2())
                .count();
        if(!Double.isFinite(frac)) {
            frac = 0.0;
        }
        this.fractionOfTopScoringContacts = frac;
    }

    public int getResidueIdentifier1() {
        return residueIdentifier1;
    }

    public String getAa1() {
        return aa1;
    }

    public String getAa2() {
        return aa2;
    }

    public int getResidueIdentifier2() {
        return residueIdentifier2;
    }

    public ContactDistanceBin getContactDistanceBin() {
        return contactDistanceBin;
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

    public double getAverageRmsdIncrease() {
        return averageRmsdIncrease;
    }

    public double getAverageTmScoreIncrease() {
        return averageTmScoreIncrease;
    }

    public double getAverageQIncrease() {
        return averageQIncrease;
    }

    public double getMaximumRmsdIncrease() {
        return maximumRmsdIncrease;
    }

    public double getMaximumTmScoreIncrease() {
        return maximumTmScoreIncrease;
    }

    public double getMaximumQIncrease() {
        return maximumQIncrease;
    }

    public boolean isEarlyFoldingResidue() {
        return isEarlyFoldingResidue;
    }

    public boolean isEarlyFoldingContact() {
        return isEarlyFoldingContact;
    }

    public double getPlmScore() {
        return plmScore;
    }

    public void setPlmScore(double plmScore) {
        this.plmScore = plmScore;
    }

    public double getAverageRmsdIncreaseZScore() {
        return averageRmsdIncreaseZScore;
    }

    public double getMaximumRmsdIncreaseZScore() {
        return maximumRmsdIncreaseZScore;
    }

    public double getFractionOfTopScoringContacts() {
        return fractionOfTopScoringContacts;
    }

    public String getCsvLine() {
        return residueIdentifier1 + "," +
                residueIdentifier2 + "," +
                contactDistanceBin + "," +
                StandardFormat.format(baselineRmsd) + "," +
                StandardFormat.format(baselineTmScore) + "," +
                StandardFormat.format(baselineQ) + "," +
                StandardFormat.format(averageRmsdIncrease) + "," +
                StandardFormat.format(averageTmScoreIncrease) + "," +
                StandardFormat.format(averageQIncrease) + "," +
                StandardFormat.format(maximumRmsdIncrease) + "," +
                StandardFormat.format(maximumTmScoreIncrease) + "," +
                StandardFormat.format(maximumQIncrease) + "," +
                isEarlyFoldingResidue + "," +
                isEarlyFoldingContact;
    }
}
