package pmw;

import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

/**
 * Created by S on 24.10.2016.
 */
public class SequenceMotifAnnotatorFunctionalTest {
    private ProteinParser proteinParser;
    private Protein protein1acj;
    private Protein protein1brr;

    @Before
    public void setup() throws IOException {
        proteinParser = new ProteinParser();
        protein1acj = proteinParser.parseProteinById("1acj");
        protein1brr = proteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldAnnotateSequenceMotifs() {
        SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();
        sequenceMotifAnnotator.process(protein1acj);
        //TODO can map entries be null?
        protein1acj.residues().map(residue -> residue.getFeature(List.class,
                SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF.name())).forEach(System.out::println);
    }
}
