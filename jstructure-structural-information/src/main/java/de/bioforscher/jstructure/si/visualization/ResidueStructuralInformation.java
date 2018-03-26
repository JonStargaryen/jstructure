package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

public class ResidueStructuralInformation {
    private final ResidueIdentifier residueIdentifier;
    private final double sumAverageRmsdIncrease;
    private final double sumAverageTmScoreIncrease;
    private final double sumAverageQIncrease;
    private final double sumMaxRmsdIncrease;
    private final double sumMaxTmScoreIncrease;
    private final double sumMaxQIncrease;
    private final boolean isEarlyFoldingResidue;

    public ResidueStructuralInformation(ResidueIdentifier residueIdentifier,
                                        double sumAverageRmsdIncrease,
                                        double sumAverageTmScoreIncrease,
                                        double sumAverageQIncrease,
                                        double sumMaxRmsdIncrease,
                                        double sumMaxTmScoreIncrease,
                                        double sumMaxQIncrease,
                                        boolean isEarlyFoldingResidue) {
        this.residueIdentifier = residueIdentifier;
        this.sumAverageRmsdIncrease = sumAverageRmsdIncrease;
        this.sumAverageTmScoreIncrease = sumAverageTmScoreIncrease;
        this.sumAverageQIncrease = sumAverageQIncrease;
        this.sumMaxRmsdIncrease = sumMaxRmsdIncrease;
        this.sumMaxTmScoreIncrease = sumMaxTmScoreIncrease;
        this.sumMaxQIncrease = sumMaxQIncrease;
        this.isEarlyFoldingResidue = isEarlyFoldingResidue;
    }

    public ResidueIdentifier getResidueIdentifier() {
        return residueIdentifier;
    }

    public double getSumAverageRmsdIncrease() {
        return sumAverageRmsdIncrease;
    }

    public double getSumAverageTmScoreIncrease() {
        return sumAverageTmScoreIncrease;
    }

    public double getSumAverageQIncrease() {
        return sumAverageQIncrease;
    }

    public double getSumMaxRmsdIncrease() {
        return sumMaxRmsdIncrease;
    }

    public double getSumMaxTmScoreIncrease() {
        return sumMaxTmScoreIncrease;
    }

    public double getSumMaxQIncrease() {
        return sumMaxQIncrease;
    }

    public boolean isEarlyFoldingResidue() {
        return isEarlyFoldingResidue;
    }

    public String getCsvLine() {
        return residueIdentifier + "," +
                StandardFormat.format(sumAverageRmsdIncrease) + "," +
                StandardFormat.format(sumAverageTmScoreIncrease) + "," +
                StandardFormat.format(sumAverageQIncrease) + "," +
                StandardFormat.format(sumMaxRmsdIncrease) + "," +
                StandardFormat.format(sumMaxTmScoreIncrease) + "," +
                StandardFormat.format(sumMaxQIncrease) + "," +
                isEarlyFoldingResidue;
    }
}
