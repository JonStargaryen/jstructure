package design;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import org.junit.Test;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by S on 24.11.2016.
 */
public class SequenceBoxRendererTest {
    @Test
    public void shouldRenderSequenceBox() throws IOException {
        //TODO make paths point to test/resources
        List<String> sequences = Files.readAllLines(Paths.get("d:/GG4-tm.seq"));

        SequenceBoxRenderer.renderSequenceBox(sequences, Paths.get("d:/gg4-tm.html"), SequenceBoxRenderer.IDENTITY_MAPPING);
        SequenceBoxRenderer.renderSequenceBox(sequences, Paths.get("d:/gg4-tm-gut.html"), SequenceBoxRenderer.GUTTERIDGE_MAPPING);
    }

    @Test
    public void shouldComposeHtmlOfAllStructures() throws IOException {
        List<String> topologies = Arrays.asList("tm", "ntm");

        String output = Stream.of(SequenceMotifDefinition.values())
                .map(sequenceMotifDefinition -> topologies.stream()
                        .map(topology -> {
                            try {
                                List<String> sequences = Files.readAllLines(Paths.get(DesignConstants.EXTRACTED_SEQUENCES_BY_TOPOLOGY_DIR + topology + "/" + sequenceMotifDefinition + DesignConstants.SEQUENCE_SUFFIX));
                                String rawOutput = SequenceBoxRenderer.renderSequenceBox(sequences, null, SequenceBoxRenderer.IDENTITY_MAPPING);
                                return "<h1>" + sequenceMotifDefinition + " - " + topology + "</h1>" + rawOutput;
                            } catch (IOException e) {
                                throw new UncheckedIOException(e);
                            }
                        })
                        .collect(Collectors.joining("</div><div class=\"topology-container\">",
                                "<div class=\"topology-container\">",
                                "</div>")))
                .collect(Collectors.joining("</div><div class=\"motif-container\">",
                        "<div class=\"motif-container\">",
                        "</div>"));

        Files.write(Paths.get("d:/output.html"), output.getBytes());
    }
}