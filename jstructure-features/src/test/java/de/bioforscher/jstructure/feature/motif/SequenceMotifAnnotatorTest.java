package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;

/**
 * The sequence motif annotator test.
 * Created by S on 24.10.2016.
 */
public class SequenceMotifAnnotatorTest {
    private Structure protein;

    @Before
    public void setup() throws IOException {
        protein = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
    }

    @Test
    public void shouldAnnotateSequenceMotifs() {
        SequenceMotifAnnotator sequenceMotifAnnotator = new SequenceMotifAnnotator();
        sequenceMotifAnnotator.process(protein);
        protein.aminoAcids()
                .forEach(aminoAcid -> Assert.assertNotNull(aminoAcid.getFeature(SequenceMotifContainer.class)));
        SequenceMotifContainer container = protein.getFeature(SequenceMotifContainer.class);
        Assert.assertEquals("GG4 count does not match", 5, count(container, SequenceMotifDefinition.GG4));
        Assert.assertEquals("AL6 count does not match", 2, count(container, SequenceMotifDefinition.AL6));
        Assert.assertEquals("IL4 count does not match", 1, count(container, SequenceMotifDefinition.IL4));
    }

    private long count(SequenceMotifContainer sequenceMotifContainer, SequenceMotifDefinition sequenceMotifDefinition) {
        return sequenceMotifContainer.getSequenceMotifs()
                .stream()
                .filter(sequenceMotif -> sequenceMotif.getMotifDefinition() == sequenceMotifDefinition)
                .count();
    }
}
