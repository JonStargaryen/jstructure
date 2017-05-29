package de.bioforscher.jstructure.model.structure;

import org.junit.Test;

/**
 * Test the group prototype parser.
 * Created by S on 24.05.2017.
 */
public class GroupPrototypeParserTest {
    @Test
    public void shouldParseStandardAminoAcid() {
        System.out.println(GroupPrototypeParser.getInstance().getPrototype("ALA"));
    }
}