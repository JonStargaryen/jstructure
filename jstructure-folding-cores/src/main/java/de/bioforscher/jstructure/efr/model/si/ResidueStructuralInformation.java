package de.bioforscher.jstructure.efr.model.si;

import de.bioforscher.jstructure.StandardFormat;

import java.util.List;

public class ResidueStructuralInformation {
    private final int residueIdentifier;
    private final String aa;
    private final double averageRmsdIncrease;
    private final double averageTmScoreIncrease;
    private final double averageQIncrease;
    private final double maximumRmsdIncrease;
    private final double maximumTmScoreIncrease;
    private final double maximumQIncrease;
    private final boolean isEarlyFoldingResidue;
    private final double eccount;
    private final double cumstrength;
    private final double conservation;
    private final double averageRmsdIncreaseZScore;
    private final double maximumRmsdIncreaseZScore;
    private final double fractionOfTopScoringContacts;
    private final double sumRmsdIncrease;
    private final List<Double> rawValues;

    public ResidueStructuralInformation(int residueIdentifier,
                                        String aa,
                                        double averageRmsdIncrease,
                                        double averageTmScoreIncrease,
                                        double averageQIncrease,
                                        double maximumRmsdIncrease,
                                        double maximumTmScoreIncrease,
                                        double maximumQIncrease,
                                        boolean isEarlyFoldingResidue,
                                        double eccount,
                                        double cumstrength,
                                        double conservation,
                                        double averageRmsdIncreaseZScore,
                                        double maximumRmsdIncreaseZScore,
                                        double fractionOfTopScoringContacts,
                                        double sumRmsdIncrease,
                                        List<Double> rawValues) {
        this.residueIdentifier = residueIdentifier;
        this.aa = aa;
        this.averageRmsdIncrease = averageRmsdIncrease;
        this.averageTmScoreIncrease = averageTmScoreIncrease;
        this.averageQIncrease = averageQIncrease;
        this.maximumRmsdIncrease = maximumRmsdIncrease;
        this.maximumTmScoreIncrease = maximumTmScoreIncrease;
        this.maximumQIncrease = maximumQIncrease;
        this.isEarlyFoldingResidue = isEarlyFoldingResidue;
        this.eccount = eccount;
        this.cumstrength = cumstrength;
        this.conservation = conservation;
        this.averageRmsdIncreaseZScore = averageRmsdIncreaseZScore;
        this.maximumRmsdIncreaseZScore = maximumRmsdIncreaseZScore;
        this.fractionOfTopScoringContacts = fractionOfTopScoringContacts;
        this.sumRmsdIncrease = sumRmsdIncrease;
        this.rawValues = rawValues;
    }

    public int getResidueIdentifier() {
        return residueIdentifier;
    }

    public String getAa() {
        return aa;
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

    public double getEccount() {
        return eccount;
    }

    public double getCumstrength() {
        return cumstrength;
    }

    public double getConservation() {
        return conservation;
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

    public double getSumRmsdIncrease() {
        return sumRmsdIncrease;
    }

    public List<Double> getRawValues() {
        return rawValues;
    }

    public String getCsvLine() {
        return residueIdentifier + "," +
                StandardFormat.format(averageRmsdIncrease) + "," +
                StandardFormat.format(averageTmScoreIncrease) + "," +
                StandardFormat.format(averageQIncrease) + "," +
                StandardFormat.format(maximumRmsdIncrease) + "," +
                StandardFormat.format(maximumTmScoreIncrease) + "," +
                StandardFormat.format(maximumQIncrease) + "," +
                isEarlyFoldingResidue;
    }

    @Override
    public String toString() {
        return "ResidueStructuralInformation{" +
                "residueIdentifier=" + residueIdentifier +
                ", aa='" + aa + '\'' +
                ", averageRmsdIncrease=" + averageRmsdIncrease +
                ", averageTmScoreIncrease=" + averageTmScoreIncrease +
                ", averageQIncrease=" + averageQIncrease +
                ", maximumRmsdIncrease=" + maximumRmsdIncrease +
                ", maximumTmScoreIncrease=" + maximumTmScoreIncrease +
                ", maximumQIncrease=" + maximumQIncrease +
                ", isEarlyFoldingResidue=" + isEarlyFoldingResidue +
                ", eccount=" + eccount +
                ", cumstrength=" + cumstrength +
                ", conservation=" + conservation +
                ", averageRmsdIncreaseZScore=" + averageRmsdIncreaseZScore +
                ", maximumRmsdIncreaseZScore=" + maximumRmsdIncreaseZScore +
                ", fractionOfTopScoringContacts=" + fractionOfTopScoringContacts +
                ", sumRmsdIncrease=" + sumRmsdIncrease +
                '}';
    }
}
