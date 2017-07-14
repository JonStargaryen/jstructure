package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Map;

/**
 * Describes the degree of sequential conservation of an amino acid.
 * Created by bittrich on 7/13/17.
 */
public class SequenceConservationProfile extends FeatureContainerEntry {
    private final Map<AminoAcid.Family, Double> frequenceTable;
    private final double deletionFrequency;

    public SequenceConservationProfile(Map<AminoAcid.Family, Double> frequencyTable,
                                       double deletionFrequency) {
        super(null);
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
