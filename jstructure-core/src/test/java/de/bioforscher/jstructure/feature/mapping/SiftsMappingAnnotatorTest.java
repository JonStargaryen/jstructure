package de.bioforscher.jstructure.feature.mapping;

import de.bioforscher.jstructure.model.identifier.PdbChainId;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

/**
 * Test the SIFTS-mapping.
 * Created by bittrich on 4/6/17.
 */
public class SiftsMappingAnnotatorTest {
    private SiftsMappingAnnotator annotator = new SiftsMappingAnnotator();

    @Test
    public void shouldHandle5tgz() {
        Protein protein = ProteinParser.source("5tgz").parse();
        annotator.process(protein);
    }

    @Test
    public void shouldDetermineUniProtNumbering() {
        Protein protein = ProteinParser.source("1acj").parse();
        annotator.process(protein);

        ChainMapping chainSiftsMapping = protein.getChains().get(0).getFeatureContainer().getFeature(ChainMapping.class);
        Assert.assertEquals("P04058", chainSiftsMapping.getUniProtId());
        Assert.assertEquals("3.1.1.7", chainSiftsMapping.getEcNumber());
        Assert.assertEquals("PF00135", chainSiftsMapping.getPfamId());

        protein.aminoAcids().forEach(group -> Assert.assertTrue(group.getFeatureContainer().getFeatureOptional(ResidueMapping.class).isPresent()));
    }

    @Test
    public void testXmlDownload() {
        Document document = annotator.downloadXml("1acj");
        Assert.assertTrue(document.html().startsWith("<?xml"));
    }

    @Test
    public void testMappingToParentElement() {
        Document document = annotator.downloadXml("1acj");
        Element element = annotator.mapToDescribingElement(document, "A", "4");
        Assert.assertTrue(element.html().contains("P04058"));
        Assert.assertTrue(element.html().contains("dbResNum=\"25\""));
    }

    @Test
    public void testMappingToUniProt() {
        Document document = annotator.downloadXml("1acj");
        ResidueMapping uniProtResidueNumber = annotator.mapGroup(document, "A", 4);
        Assert.assertEquals("25", uniProtResidueNumber.getUniProtResidueNumber());
    }
}