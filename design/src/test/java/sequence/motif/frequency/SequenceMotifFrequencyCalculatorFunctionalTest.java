package sequence.motif.frequency;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;
import design.sequence.motif.frequency.SequenceMotifFrequencyCalculator;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by S on 24.11.2016.
 */
public class SequenceMotifFrequencyCalculatorFunctionalTest {
    @Test
    public void shouldComposeHtmlOfAllStructures() throws IOException {
        final String basePath = DesignConstants.EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR;
        final String ext = DesignConstants.SEQUENCE_SUFFIX;

        Map<SequenceMotifDefinition, List<SequenceMotifFrequencyCalculator>> frequencies =
                Stream.of(SequenceMotifDefinition.values())
                .peek(System.out::println)
                .map(motif -> new SequenceMotifFrequencyCalculator(motif, Paths.get(basePath + "tm/" + motif + ext),
                        Paths.get(basePath + "ntm/" + motif + ext), Paths.get(basePath + "trans/" + motif + ext)))
                .collect(Collectors.groupingBy(SequenceMotifFrequencyCalculator::getSequenceMotif));

        ObjectMapper objectMapper = new ObjectMapper();
        System.out.println(objectMapper
                // indent output
                .writerWithDefaultPrettyPrinter()
                .writeValueAsString(frequencies)
                // trim double precision - regex 101:
                // ?<= - this group has to be present, but will not be replaced
                .replaceAll("(?<=\\d\\.\\d{4})\\d*", "")
        );
    }
}
