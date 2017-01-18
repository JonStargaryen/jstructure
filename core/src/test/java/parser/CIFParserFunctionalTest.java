package parser;

import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.CIFParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * The functional test for the CIF parser.
 * Created by bittrich on 1/18/17.
 */
public class CIFParserFunctionalTest {
    @Test
    public void shouldAnnotateMPH() {
        GroupInformation groupInformation = CIFParser.parseLigandInformation("MPH");
        //TODO this should fail probably
        Assert.assertFalse(groupInformation.isAminoAcid());
    }
}
