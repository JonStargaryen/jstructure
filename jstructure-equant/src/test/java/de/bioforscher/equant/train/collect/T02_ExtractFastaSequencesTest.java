package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

public class T02_ExtractFastaSequencesTest {
    @Test
    @Ignore
    public void shouldReportAminoAcidSequence() {
        Structure structure = StructureParser.source(EquantConstants.CASP9_DIRECTORY.resolve("targets").resolve("T0643.pdb")).parse();
        Assert.assertTrue(structure.getAminoAcidSequence().length() > 0);
    }
}