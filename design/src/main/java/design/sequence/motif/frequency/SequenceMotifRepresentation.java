package design.sequence.motif.frequency;

import de.bioforscher.jstructure.model.structure.AminoAcid;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Counts the occurrences of individual amino acids in an ensemble of protein sequences.
 * Created by S on 13.12.2016.
 */
public class SequenceMotifRepresentation {
    // list: sequence motif position, map: aa-grouping -> frequency
    private final List<Map<String, Double>> frequencies;

    public SequenceMotifRepresentation(Path sequenceFile) throws IOException {
        this(Files.readAllLines(sequenceFile));
    }

    public SequenceMotifRepresentation(List<String> sequences) {
        this.frequencies = countFrequencies(sequences);
    }

    /**
     * Computes the delta of frequencies between tm and ntm regions.
     *
     * @param transmembrane    the tm frequencies
     * @param nonTransmembrane the ntm frequencies
     */
    public SequenceMotifRepresentation(SequenceMotifRepresentation transmembrane, SequenceMotifRepresentation nonTransmembrane) {
        this.frequencies = new ArrayList<>();

        // merge key sets, so the difference can and will include all initial entries
        Set<String> commonKeySet = Stream.of(transmembrane.frequencies, nonTransmembrane.frequencies)
                .flatMap(Collection::stream)
                .map(Map::keySet)
                .flatMap(Set::stream)
                .collect(Collectors.toSet());

        for (int position = 0; position < transmembrane.frequencies.size(); position++) {
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

    private List<Map<String, Double>> countFrequencies(List<String> sequences) {
        List<Map<String, Double>> aminoAcidFrequencies = new ArrayList<>();
        //TODO when did this '-1' happen?
        final int sequenceLength = sequences.get(0).length() - 1;
        final int sequenceCount = sequences.size();

        // traverse all variable position of sequence motif and count relative occurrence of each position by the
        // specified mapping rule
        IntStream.range(1, sequenceLength - 1)
                .forEach(index -> {
                    Map<String, Double> frequency = sequences.stream()
                            .filter(sequence -> sequence.length() > index)
                            .map(sequence -> String.valueOf(sequence.charAt(index)))
                            .map(AminoAcid::valueOfIgnoreCase)
                            .collect(Collectors.groupingBy(AminoAcid::getOneLetterCode, Collectors.reducing(0.0, e -> 1.0 / sequenceCount,
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