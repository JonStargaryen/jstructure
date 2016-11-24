package design;

import de.bioforscher.jstructure.model.structure.AminoAcid;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by S on 24.11.2016.
 */
public class SequenceMotifRenderer {
    enum Topology {
        TRANSMEMBRANE,
        NON_TRANSMEMBRANE,
        TRANSITION
    }
    private final Path tmSequences;
    private final Path nTmSequences;
    private final Path transSequences;

    private final Function<AminoAcid, String> mappingRule;
    // perform only a mapping to the enums name
    public static final Function<AminoAcid, String> IDENTITY_MAPPING = AminoAcid::getOneLetterCode;
    // perform mapping to the Gutteridge grouping
    public static final Function<AminoAcid, String> GUTTERIDGE_MAPPING = aminoAcid -> aminoAcid.getGutteridgeGrouping().name();

    public SequenceMotifRenderer(Path tmSequences, Path nTmSequences, Path transSequences, Function<AminoAcid, String> mappingRule) {
        this.tmSequences = tmSequences;
        this.nTmSequences = nTmSequences;
        this.transSequences = transSequences;
        this.mappingRule = mappingRule;
    }

    public SequenceMotifRenderer(Path tmSequences, Path nTmSequences, Path transSequences) {
        this(tmSequences, nTmSequences, transSequences, IDENTITY_MAPPING);
    }

    class SequenceMotifRepresentation {
        SequenceMotifRepresentation(Path sequenceFile) {

        }
    }

    public List<Map<String, Double>> countFrequencies(List<String> sequences) {
        List<Map<String, Double>> aminoAcidFrequencies = new ArrayList<>();
        final int sequenceLength = sequences.get(0).length();
        final int sequenceCount = sequences.size();

        // traverse all variable position of sequence motif and count relative occurrence of each position by the
        // specified mapping rule
        IntStream.range(1, sequenceLength - 1).forEach(index -> aminoAcidFrequencies.add(sequences.stream()
                .filter(sequence -> sequence.length() > index)
                .map(sequence -> String.valueOf(sequence.charAt(index)))
                .map(AminoAcid::valueOfIgnoreCase)
                .collect(Collectors.groupingBy(mappingRule, Collectors.reducing(0.0, e -> 1.0 / sequenceCount,
                        Double::sum)))));

        return aminoAcidFrequencies;
    }
}
