package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * The sequence motif annotator test.
 * Created by S on 24.10.2016.
 */
public class SequenceMotifAnnotatorTest {
    private Protein protein1acj;

    @Before
    public void setup() throws IOException {
        protein1acj = ProteinParser.source("1acj").parse();
    }

    @Test
    public void shouldAnnotateSequenceMotifs() {

// TODO need sophisticated tests
        SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();
        sequenceMotifAnnotator.process(protein1acj);
//        Selection.on(protein1acj)
//                .aminoAcids()
//                .asFilteredGroups()
//                .map(residue -> residue.getFeature(List.class,
//                        SequenceMotifAnnotator.SEQUENCE_MOTIF))
//                .forEach(System.out::println);
    }
}
