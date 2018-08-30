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
    private int couplingRank;
    private final List<Double> rawValues;
    private boolean top02;
    private boolean top04;
    private boolean top06;
    private boolean top08;
    private boolean top10;
    private boolean top12;
    private boolean top14;
    private boolean top16;
    private boolean hydrogenBond;
    private boolean hydrophobicInteraction;

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
                                        List<ReconstructionStructuralInformation> topScoringReconstructions,
                                        List<Double> rawValues) {
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
        this.rawValues = rawValues;
    }

    public boolean isHydrogenBond() {
        return hydrogenBond;
    }

    public boolean isHydrophobicInteraction() {
        return hydrophobicInteraction;
    }

    public int getCouplingRank() {
        return couplingRank;
    }

    public void setCouplingRank(int couplingRank) {
        this.couplingRank = couplingRank;
    }

    public List<Double> getRawValues() {
        return rawValues;
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

    public void markAsTopScoringContact02() {
        this.top02 = true;
    }

    public void markAsTopScoringContact04() {
        this.top04 = true;
    }

    public void markAsTopScoringContact06() {
        this.top06 = true;
    }

    public void markAsTopScoringContact08() {
        this.top08 = true;
    }

    public void markAsTopScoringContact10() {
        this.top10 = true;
    }

    public void markAsTopScoringContact12() {
        this.top12 = true;
    }

    public void markAsTopScoringContact14() {
        this.top14 = true;
    }

    public void markAsTopScoringContact16() {
        this.top16 = true;
    }

    public void markAsHydrogenBond() {
        hydrogenBond = true;
    }

    public void markAsHydrophobicInteraction() {
        hydrophobicInteraction = true;
    }

    public boolean isTop02() {
        return top02;
    }

    public boolean isTop04() {
        return top04;
    }

    public boolean isTop06() {
        return top06;
    }

    public boolean isTop08() {
        return top08;
    }

    public boolean isTop10() {
        return top10;
    }

    public boolean isTop12() {
        return top12;
    }

    public boolean isTop14() {
        return top14;
    }

    public boolean isTop16() {
        return top16;
    }

    public boolean istop02() {
        return top02;
    }
}
