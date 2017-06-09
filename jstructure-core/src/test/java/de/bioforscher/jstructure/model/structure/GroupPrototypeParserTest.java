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

    @Test
    public void shouldGetDocumentForAlanine() {
        System.out.println(GroupPrototypeParser.getDocument("ALA"));
    }

    @Test
    public void shouldHandleUnknownLigand() {
        System.out.println(GroupPrototypeParser.getInstance().getPrototype("UNL"));
    }

    @Test
    public void s3TAK() {
        System.out.println("3TAK");
    }
}