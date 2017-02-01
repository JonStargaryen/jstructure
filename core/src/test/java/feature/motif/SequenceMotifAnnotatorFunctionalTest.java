package feature.motif;

import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

/**
 * The sequence motif annotator test.
 * TODO need sophisticated tests
 * Created by S on 24.10.2016.
 */
public class SequenceMotifAnnotatorFunctionalTest {
    private Protein protein1acj;
    private Protein protein1brr;

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.parseProteinById("1acj");
        protein1brr = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldAnnotateSequenceMotifs() {
        SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();
        sequenceMotifAnnotator.process(protein1acj);
        Selection.on(protein1acj)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getFeature(List.class,
                        SequenceMotifAnnotator.SEQUENCE_MOTIF))
                .forEach(System.out::println);
    }
}
