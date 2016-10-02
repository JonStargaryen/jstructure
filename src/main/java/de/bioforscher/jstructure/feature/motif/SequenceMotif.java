package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.feature.ComplexFeatureValue;
import de.bioforscher.jstructure.model.structure.Residue;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif implements ComplexFeatureValue {
    private final SequenceMotifDefinition motifDefinition;
    private final Residue startResidue;
    private final Residue endResidue;
    private final String sequence;

    public SequenceMotif(SequenceMotifDefinition candidate, Residue startResidue, Residue endResidue, String sequence) {
        this.motifDefinition = candidate;
        this.startResidue = startResidue;
        this.endResidue = endResidue;
        this.sequence = sequence;
    }

    public String getSequence() {
        return sequence;
    }

    public Residue getEndResidue() {
        return endResidue;
    }

    public Residue getStartResidue() {
        return startResidue;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " motifDefinition='" + motifDefinition + "' sequence='" + sequence +
                "' startResidue='" + startResidue.getName() + "' endResidue='" + endResidue.getName() + "'";
    }
}
