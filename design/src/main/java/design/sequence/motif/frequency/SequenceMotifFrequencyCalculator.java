package design.sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.AminoAcid;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Computes frequency statistics on a set of 3 sequence files describing a
 * {@link de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition} in 3 different topologies.
 * Created by S on 24.11.2016.
 */
public class SequenceMotifFrequencyCalculator {
    private final Function<AminoAcid, String> mappingRule;
    // perform only a mapping to the enums name
    public static final Function<AminoAcid, String> IDENTITY_MAPPING = AminoAcid::getOneLetterCode;
    // perform mapping to the Gutteridge grouping
    public static final Function<AminoAcid, String> GUTTERIDGE_MAPPING = aminoAcid -> aminoAcid.getGutteridgeGrouping().name();

    private final SequenceMotifRepresentation transmembraneFrequencies;
    private final SequenceMotifRepresentation nonTransmembraneFrequencies;
    private final SequenceMotifRepresentation transistionFrequencies;
    // difference between transmembrane and non-transmembrane frequencies: f_tm - f_ntm
    private final SequenceMotifRepresentation deltaMembraneFrequencies;
    private final SequenceMotifDefinition sequenceMotif;

    public SequenceMotifFrequencyCalculator(SequenceMotifDefinition motif, Path tmSequences, Path nTmSequences,
                                            Path transSequences, Function<AminoAcid, String> mappingRule) {
        this.mappingRule = mappingRule;
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

    public SequenceMotifFrequencyCalculator(SequenceMotifDefinition motif, Path tmSequences, Path nTmSequences, Path transSequences) {
        this(motif, tmSequences, nTmSequences, transSequences, IDENTITY_MAPPING);
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

    public class SequenceMotifRepresentation {
        // list: sequence motif position, map: aa-grouping -> frequency
        private final List<Map<String, Double>> frequencies;

        SequenceMotifRepresentation(Path sequenceFile) throws IOException {
            this(Files.readAllLines(sequenceFile));
        }

        SequenceMotifRepresentation(List<String> sequences) {
            this.frequencies = countFrequencies(sequences);
        }

        /**
         * Computes the delta of frequencies between tm and ntm regions.
         * @param transmembrane the tm frequencies
         * @param nonTransmembrane the ntm frequencies
         */
        SequenceMotifRepresentation(SequenceMotifRepresentation transmembrane, SequenceMotifRepresentation nonTransmembrane) {
            this.frequencies = new ArrayList<>();

            // merge key sets, so the difference can and will include all initial entries
            Set<String> commonKeySet = Stream.of(transmembrane.frequencies, nonTransmembrane.frequencies)
                    .flatMap(Collection::stream)
                    .map(Map::keySet)
                    .flatMap(Set::stream)
                    .collect(Collectors.toSet());

            for(int position = 0; position < transmembrane.frequencies.size(); position++) {
                Map<String, Double> tmPositionFrequencies = transmembrane.frequencies.get(position);
                Map<String, Double> ntmPositionFrequencies = nonTransmembrane.frequencies.get(position);
                Map<String, Double> positionFrequencies = new HashMap<>();
                commonKeySet.forEach(key -> positionFrequencies.put(key, tmPositionFrequencies.getOrDefault(key, 0.0)
                        - ntmPositionFrequencies.getOrDefault(key, 0.0)));
                frequencies.add(positionFrequencies);
            }
        }

        public List<Map<String, Double>> getFrequencies() {
            return frequencies;
        }
    }

    private List<Map<String, Double>> countFrequencies(List<String> sequences) {
        List<Map<String, Double>> aminoAcidFrequencies = new ArrayList<>();
        final int sequenceLength = sequences.get(0).length();
        final int sequenceCount = sequences.size();

        // traverse all variable position of sequence motif and count relative occurrence of each position by the
        // specified mapping rule
        IntStream.range(1, sequenceLength - 1)
                .forEach(index -> {
                    Map<String, Double> frequency = sequences.stream()
                            .filter(sequence -> sequence.length() > index)
                            .map(sequence -> String.valueOf(sequence.charAt(index)))
                            .map(AminoAcid::valueOfIgnoreCase)
                            .collect(Collectors.groupingBy(mappingRule, Collectors.reducing(0.0, e -> 1.0 / sequenceCount,
                                    Double::sum)));

                    //TODO move/make abstract, other mappings will fail here
                    Stream.of(AminoAcid.values())
                            .map(AminoAcid::getOneLetterCode)
                            .filter(olc -> !olc.equals(AminoAcid.UNKNOWN.getOneLetterCode()))
                            .forEach(olc -> frequency.putIfAbsent(olc, 0.0));

                    aminoAcidFrequencies.add(frequency);
                });

        return aminoAcidFrequencies;
    }
}
