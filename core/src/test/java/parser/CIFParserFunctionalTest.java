package parser;

import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.CIFParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test functionality of the CIF-parser.
 * Created by bittrich on 3/29/17.
 */
public class CIFParserFunctionalTest {
    @Test
    public void shouldParseNameOnMultipleLines() {
        GroupInformation groupInformation = CIFParser.parseLigandInformation("UMP");
        Assert.assertEquals("2'-DEOXYURIDINE 5'-MONOPHOSPHATE", groupInformation.getName());
    }
}
