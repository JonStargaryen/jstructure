package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;

public class ContactStructuralInformation {
    private final ResidueIdentifier residueIdentifier1;
    private final ResidueIdentifier residueIdentifier2;
    private final ContactDistanceBin contactDistanceBin;
    private final double averageRmsdIncrease;
    private final double averageTmScoreIncrease;
    private final double averageQIncrease;
    private final double maximumRmsdIncrease;
    private final double maximumTmScoreIncrease;
    private final double maximumQIncrease;
    private final boolean isEarlyFoldingResidue;
    private final boolean isEarlyFoldingContact;

    public ContactStructuralInformation(ResidueIdentifier residueIdentifier1,
                                        ResidueIdentifier residueIdentifier2,
                                        ContactDistanceBin contactDistanceBin,
                                        double averageRmsdIncrease,
                                        double averageTmScoreIncrease,
                                        double averageQIncrease,
                                        double maximumRmsdIncrease,
                                        double maximumTmScoreIncrease,
                                        double maximumQIncrease,
                                        boolean isEarlyFoldingResidue,
                                        boolean isEarlyFoldingContact) {
        this.residueIdentifier1 = residueIdentifier1;
        this.residueIdentifier2 = residueIdentifier2;
        this.contactDistanceBin = contactDistanceBin;
        this.averageRmsdIncrease = averageRmsdIncrease;
        this.averageTmScoreIncrease = averageTmScoreIncrease;
        this.averageQIncrease = averageQIncrease;
        this.maximumRmsdIncrease = maximumRmsdIncrease;
        this.maximumTmScoreIncrease = maximumTmScoreIncrease;
        this.maximumQIncrease = maximumQIncrease;
        this.isEarlyFoldingResidue = isEarlyFoldingResidue;
        this.isEarlyFoldingContact = isEarlyFoldingContact;
    }

    public ResidueIdentifier getResidueIdentifier1() {
        return residueIdentifier1;
    }

    public ResidueIdentifier getResidueIdentifier2() {
        return residueIdentifier2;
    }

    public ContactDistanceBin getContactDistanceBin() {
        return contactDistanceBin;
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

    public String getCsvLine() {
        return residueIdentifier1 + "," +
                residueIdentifier2 + "," +
                contactDistanceBin + "," +
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
