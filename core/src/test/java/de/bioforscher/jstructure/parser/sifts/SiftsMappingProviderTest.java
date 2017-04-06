package de.bioforscher.jstructure.parser.sifts;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the SIFTS-mapping.
 * Created by bittrich on 4/6/17.
 */
public class SiftsMappingProviderTest {
    private SiftsMappingProvider siftsMappingProvider;

    @Before
    public void setup() {
        siftsMappingProvider = new SiftsMappingProvider();
    }

    @Test
    public void shouldDetermineUniProtNumbering() {
        Protein protein = ProteinParser.source("1acj").parse();
        siftsMappingProvider.process(protein);

        ChainSiftsMapping chainSiftsMapping = protein.getChains().get(0).getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING);
        Assert.assertEquals("P04058", chainSiftsMapping.getUniProtId());
        Assert.assertEquals("3.1.1.7", chainSiftsMapping.getEcNumber());
        Assert.assertEquals("PF00135", chainSiftsMapping.getPfam());

        protein.aminoAcids().forEach(group -> Assert.assertTrue(group.getFeatureMap().containsKey(SiftsMappingProvider.SIFTS_MAPPING)));
    }

    @Test
    public void testXmlDownload() {
        Document document = siftsMappingProvider.downloadXml("1acj");
        Assert.assertTrue(document.html().startsWith("<?xml"));
    }

    @Test
    public void testMappingToParentElement() {
        Document document = siftsMappingProvider.downloadXml("1acj");
        Element element = siftsMappingProvider.mapToDescribingElement(document, "A", "4");
        Assert.assertTrue(element.html().contains("P04058"));
        Assert.assertTrue(element.html().contains("dbResNum=\"25\""));
    }

    @Test
    public void testMappingToUniProt() {
        Document document = siftsMappingProvider.downloadXml("1acj");
        ResidueSiftsMapping uniProtResidueNumber = siftsMappingProvider.mapGroup(document, "A", 4);
        Assert.assertEquals("P04058", uniProtResidueNumber.getUniProtId());
        Assert.assertEquals(25, uniProtResidueNumber.getUniProtResidueNumber());
    }
}