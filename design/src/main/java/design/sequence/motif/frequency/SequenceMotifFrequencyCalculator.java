package design.sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Path;

/**
 * Computes frequency statistics on a set of 3 sequence files describing a
 * {@link de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition} in 3 different topologies.
 * Created by S on 24.11.2016.
 */
public class SequenceMotifFrequencyCalculator {
    private final SequenceMotifRepresentation transmembraneFrequencies;
    private final SequenceMotifRepresentation nonTransmembraneFrequencies;
    private final SequenceMotifRepresentation transistionFrequencies;
    // difference between transmembrane and non-transmembrane frequencies: f_tm - f_ntm
    private final SequenceMotifRepresentation deltaMembraneFrequencies;
    private final SequenceMotifDefinition sequenceMotif;

    public SequenceMotifFrequencyCalculator(SequenceMotifDefinition motif, Path tmSequences, Path nTmSequences,
                                            Path transSequences) {
        this.sequenceMotif = motif;

        try {
            this.transmembraneFrequencies = new SequenceMotifRepresentation(tmSequences);
            this.nonTransmembraneFrequencies = new SequenceMotifRepresentation(nTmSequences);
            this.transistionFrequencies = new SequenceMotifRepresentation(transSequences);
            this.deltaMembraneFrequencies = new SequenceMotifRepresentation(transmembraneFrequencies, nonTransmembraneFrequencies);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public SequenceMotifDefinition getSequenceMotif() {
        return sequenceMotif;
    }

    public SequenceMotifRepresentation getTransmembraneFrequencies() {
        return transmembraneFrequencies;
    }

    public SequenceMotifRepresentation getNonTransmembraneFrequencies() {
        return nonTransmembraneFrequencies;
    }

    public SequenceMotifRepresentation getTransistionFrequencies() {
        return transistionFrequencies;
    }

    public SequenceMotifRepresentation getDeltaMembraneFrequencies() {
        return deltaMembraneFrequencies;
    }
}
