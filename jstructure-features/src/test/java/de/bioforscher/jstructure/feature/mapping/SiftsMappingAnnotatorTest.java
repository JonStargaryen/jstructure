package de.bioforscher.jstructure.feature.mapping;

import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test the SIFTS-mapping.
 * Created by bittrich on 4/6/17.
 */
public class SiftsMappingAnnotatorTest {
    private SiftsMappingAnnotator mappingAnnotator = new SiftsMappingAnnotator();

    @Test
    public void shouldMap4r3z_A() {
        Protein protein = ProteinParser.source("4r3z")
                .minimalParsing(true)
                .parse();
        mappingAnnotator.process(protein);
        ChainMapping chainMapping = protein.select()
                .chainName("A")
                .asChain()
                .getFeatureContainer()
                .getFeature(ChainMapping.class);
        Assert.assertEquals("PF01588", chainMapping.getPfamId());
        Assert.assertEquals("?", chainMapping.getEcNumber());
        Assert.assertEquals("Q12904", chainMapping.getUniProtId());
    }

    @Test
    public void shouldMap1o0c_A() {
        Protein protein = ProteinParser.source("1o0c")
                .minimalParsing(true)
                .parse();
        mappingAnnotator.process(protein);
        ChainMapping chainMapping = protein.select()
                .chainName("A")
                .asChain()
                .getFeatureContainer()
                .getFeature(ChainMapping.class);
        Assert.assertEquals("?", chainMapping.getPfamId());
        Assert.assertEquals("6.1.1.18", chainMapping.getEcNumber());
        Assert.assertEquals("P00962", chainMapping.getUniProtId());
    }

    @Test(expected = ComputationException.class)
    public void shouldFailForUnknownPdbId() {
        Protein protein = ProteinParser.source(TestUtils.getResourceAsInputStream("uniprot/1bs2_A_renum.pdb"))
                .parse();
        mappingAnnotator.process(protein);
    }

    @Test
    public void shouldHandle5tgl() {
        Protein protein = ProteinParser.source("5tgl").parse();
        mappingAnnotator.process(protein);
    }

    @Test
    public void shouldDetermineUniProtNumbering() {
        Protein protein = ProteinParser.source("1acj").parse();
        mappingAnnotator.process(protein);

        ChainMapping chainSiftsMapping = protein.getChains().get(0).getFeatureContainer().getFeature(ChainMapping.class);
        Assert.assertEquals("P04058", chainSiftsMapping.getUniProtId());
        Assert.assertEquals("3.1.1.7", chainSiftsMapping.getEcNumber());
        Assert.assertEquals("PF00135", chainSiftsMapping.getPfamId());

        protein.aminoAcids().forEach(group -> Assert.assertTrue(group.getFeatureContainer().getFeatureOptional(ResidueMapping.class).isPresent()));
    }

    @Test
    public void testXmlDownload() {
        Document document = mappingAnnotator.downloadXml("1acj");
        Assert.assertTrue(document.html().startsWith("<?xml"));
    }

    @Test
    public void testMappingToParentElement() {
        Document document = mappingAnnotator.downloadXml("1acj");
        Element element = mappingAnnotator.mapToDescribingElement(document, "A", "4");
        Assert.assertTrue(element.html().contains("P04058"));
        Assert.assertTrue(element.html().contains("dbResNum=\"25\""));
    }

    @Test
    public void testMappingToUniProt() {
        Document document = mappingAnnotator.downloadXml("1acj");
        ResidueMapping uniProtResidueNumber = mappingAnnotator.mapGroup(document, "A", 4);
        Assert.assertEquals("25", uniProtResidueNumber.getUniProtResidueNumber());
    }
}