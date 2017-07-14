package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Map;

/**
 * The sequence conservation profile for a given amino acid.
 * Created by bittrich on 7/12/17.
 */
@Deprecated
public class SequenceConservationProfile extends FeatureContainerEntry {
    private final Map<AminoAcid.Family, Double> frequenceTable;
    private final double deletionFrequency;

    public SequenceConservationProfile(AbstractFeatureProvider featureProvider,
                                       Map<AminoAcid.Family, Double> frequencyTable,
                                       double deletionFrequency) {
        super(featureProvider);
        this.frequenceTable = frequencyTable;
        this.deletionFrequency = deletionFrequency;
    }

    public Map<AminoAcid.Family, Double> getFrequenceTable() {
        return frequenceTable;
    }

    public double getFrequency(AminoAcid.Family aminoAcid) {
        return frequenceTable.get(aminoAcid);
    }

    public double getDeletionFrequency() {
        return deletionFrequency;
    }

    @Override
    public String toString() {
        return frequenceTable + " deletionFrequency=" + deletionFrequency;
    }
}
