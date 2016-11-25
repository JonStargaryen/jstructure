package design.sequence.motif.frequency;

import de.bioforscher.jstructure.model.structure.AminoAcid;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Converts an 'alignment' of sequences to an HTML representation.
 * Created by S on 24.11.2016.
 */
@Deprecated
public class SequenceBoxRenderer {
    private static final String SEQUENCE_BOX_TEMPLATE = "design/sequenceBoxTemplate.html";
    private String template;
    private final List<String> sequences;
    private final int sequenceLength;
    private final int sequenceCount;
    private List<Map<String, Double>> aminoAcidFrequencies;
    private final Function<AminoAcid, String> mappingRule;
    // perform only a mapping to the enums name
    public static final Function<AminoAcid, String> IDENTITY_MAPPING = AminoAcid::getOneLetterCode;
    // perform mapping to the Gutteridge grouping
    public static final Function<AminoAcid, String> GUTTERIDGE_MAPPING = aminoAcid -> aminoAcid.getGutteridgeGrouping().name();

    public SequenceBoxRenderer(List<String> sequences) {
        this(sequences, IDENTITY_MAPPING);
    }

    public SequenceBoxRenderer(List<String> sequences, Function<AminoAcid, String> mappingRule) {
        this.template = loadTemplate();
        this.sequences = sequences;
        this.sequenceCount = sequences.size();
        this.sequenceLength = sequences.get(0).length();
        this.aminoAcidFrequencies = new ArrayList<>();
        this.mappingRule = mappingRule;
    }

    public static String renderSequenceBox(List<String> sequences, Path outputHtml, Function<AminoAcid, String> mapping) throws IOException {
        SequenceBoxRenderer sequenceBoxRenderer = new SequenceBoxRenderer(sequences, mapping);
        sequenceBoxRenderer.renderSequenceBox();
        if(outputHtml != null) {
            Files.write(outputHtml, sequenceBoxRenderer.template.getBytes());
        }
        return sequenceBoxRenderer.template;
    }

    private void renderSequenceBox() {
        countOccurrences();

//        // should extract frequencies
//        aminoAcidFrequencies.stream()
//                .map(Object::toString)
//                .forEach(System.out::println);
//
//        // frequencies should sum up to 1.0
//        aminoAcidFrequencies.stream()
//                .map(map -> map.values().stream().mapToDouble(Double::valueOf).sum())
//                .forEach(System.out::println);

        composeHtml();
    }

    private void composeHtml() {
        //TODO move to mature HTML generation/processing
        template = template.replace("{{title}}", sequences.get(0))
                           .replace("{{frequencies}}", composeFrequencyTable())
        // show sequences
//                           .replace("{{sequences}}", sequences.stream()
//                                   .collect(Collectors.joining("</div>" + System.lineSeparator() + "<div>",
//                                           "<div>",
//                                           "</div>")))
        // hide sequences
                            .replace("{{sequences}}", "")
        ;
    }

    private static final DecimalFormat numberFormat = new DecimalFormat("00.000", DecimalFormatSymbols.getInstance(Locale.US));

    private String composeFrequencyTable() {
        Set<String> uniqueKeys;
        if(mappingRule.equals(IDENTITY_MAPPING) || mappingRule.equals(GUTTERIDGE_MAPPING)) {
            uniqueKeys = Stream.of(AminoAcid.values())
                    .filter(aminoAcid -> !aminoAcid.equals(AminoAcid.UNKNOWN))
                    .map(mappingRule)
                    .collect(Collectors.toSet());
        } else {
            // no standard mapping used
            uniqueKeys = aminoAcidFrequencies.stream()
                    .flatMap(map -> map.keySet().stream())
                    .collect(Collectors.toSet());
        }

        return aminoAcidFrequencies.stream()
                .map(map -> uniqueKeys.stream()
                        .map(map::get)
                        .map(value -> numberFormat.format((value != null ? value : new Double(0.0)) * 100))
                        .map(value -> "<td style=\"background-color: hsl(240, 100%, " + (100 - Double.valueOf(value)) + "%)\">" +
//                                value +
                                "</td>")
                        .collect(Collectors.joining(System.lineSeparator())))
                .collect(Collectors.joining("</tr>" + System.lineSeparator() + "<tr>",
                        "<table><tr>" +
                                uniqueKeys.stream()
                                        .map(key -> "<td class=\"header\">" + key + "</td>")
                                        .collect(Collectors.joining(System.lineSeparator()))
                                + "</tr>" + System.lineSeparator() + "<tr>",
                        "</tr></table>"));
    }

    private void countOccurrences() {
        // traverse all variable position of sequence motif and count relative occurrence of each position by the
        // specified mapping rule
        IntStream.range(1, sequenceLength - 1).forEach(index -> aminoAcidFrequencies.add(sequences.stream()
                .filter(sequence -> sequence.length() > index)
                .map(sequence -> String.valueOf(sequence.charAt(index)))
                .map(AminoAcid::valueOfIgnoreCase)
                .collect(Collectors.groupingBy(mappingRule, Collectors.reducing(0.0, e -> 1.0 / sequenceCount,
                        Double::sum)))));
    }

    private String loadTemplate() {
        try {
            return Files.lines(getResourceAsPath(SEQUENCE_BOX_TEMPLATE))
                    .collect(Collectors.joining(System.lineSeparator()));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static Path getResourceAsPath(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        Objects.requireNonNull(ccl);
        URL resource = ccl.getResource(filename);
        Objects.requireNonNull(resource);
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(resource.getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
