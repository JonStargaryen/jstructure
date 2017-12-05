package de.bioforscher.jstructure.model.structure;

import org.junit.Assert;
import org.junit.Test;

public class ElementTest {
    @Test
    public void resolveFullAtomName() throws Exception {
        Assert.assertEquals(Element.C, Element.resolveFullAtomName("CA", false));
        Assert.assertEquals(Element.Ca, Element.resolveFullAtomName("CA", true));
    }
}