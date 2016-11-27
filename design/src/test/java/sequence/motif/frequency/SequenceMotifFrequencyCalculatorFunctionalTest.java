package sequence.motif.frequency;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;
import design.sequence.motif.frequency.SequenceMotifFrequencyCalculator;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collection;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by S on 24.11.2016.
 */
public class SequenceMotifFrequencyCalculatorFunctionalTest {
    /**
     * Create container html for heatmap visualization.
     */
    @Test
    public void shouldComposeFrequencyHeatmapHTML() throws IOException {
        Map<SequenceMotifDefinition, SequenceMotifFrequencyCalculator> frequencies = loadFrequencies();

        List<String> containerIds = Stream.of(SequenceMotifDefinition.values())
                .map(Enum::name)
                .flatMap(motif -> Stream.of(motif + "-tm", motif + "-ntm", motif + "-trans", motif + "-delta"))
                .collect(Collectors.toList());

        String html = containerIds.stream()
                .map(id -> "<div class=\"" + (id.endsWith("delta") ? "delta" : "") + "\">" +
                            "<h3>" + id + "</h3>" +
                            "<div class=\"container\" id=\"" + id + "\"></div>" +
                        "</div>")
                .collect(Collectors.joining(System.lineSeparator()));

        // print html creating the container divs
        System.out.println("container html:");
        System.out.println(html);
        System.out.println();

        String data = containerIds.stream()
                // remove - for js-function & fill in data
                .map(id -> "var f" + id.replace("-", "") + " = calendarHeatmap().data(" + composeArray(id, frequencies) + ").selector('#" + id + "').colorRange(['#f4f7f7', '#79a8a9']); f" + id.replace("-", "") + "();")
                .collect(Collectors.joining(System.lineSeparator(),
                        "<script type=\"text/javascript\">",
                        System.lineSeparator() + "</script>"));

        // print data and js-code
        System.out.println("data/js");
        System.out.println(data);
    }

    private static final DecimalFormat decimalFormat = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.US));

    private String composeArray(String id, Map<SequenceMotifDefinition, SequenceMotifFrequencyCalculator> frequencies) {
        String motifName = id.split("-")[0];
        String topology = id.split("-")[1];
        SequenceMotifDefinition definition = SequenceMotifDefinition.valueOf(motifName);
        SequenceMotifFrequencyCalculator calculator = frequencies.get(definition);

        SequenceMotifFrequencyCalculator.SequenceMotifRepresentation representation;
        switch (topology) {
            case "tm":
                representation = calculator.getTransmembraneFrequencies();
                break;
            case "ntm":
                representation = calculator.getNonTransmembraneFrequencies();
                break;
            case "trans":
                representation = calculator.getTransistionFrequencies();
                break;
            case "delta":
                representation = calculator.getDeltaMembraneFrequencies();
                break;
            default:
                throw new IllegalArgumentException("unknown case: " + topology);
        }

        return representation.getFrequencies().stream()
                .map(Map::entrySet)
                .flatMap(Collection::stream)
                .map(Map.Entry::getValue)
                .map(decimalFormat::format)
                .collect(Collectors.joining(",",
                        "[",
                        "]"));
    }

    private Map<SequenceMotifDefinition, SequenceMotifFrequencyCalculator> loadFrequencies() throws IOException {
        final String basePath = DesignConstants.EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR;
        final String ext = DesignConstants.SEQUENCE_SUFFIX;

        return Stream.of(SequenceMotifDefinition.values())
                .map(motif -> new SequenceMotifFrequencyCalculator(motif, Paths.get(basePath + "tm/" + motif + ext),
                        Paths.get(basePath + "ntm/" + motif + ext), Paths.get(basePath + "trans/" + motif + ext)))
                .collect(Collectors.toMap(SequenceMotifFrequencyCalculator::getSequenceMotif, Function.identity()));
    }
}
