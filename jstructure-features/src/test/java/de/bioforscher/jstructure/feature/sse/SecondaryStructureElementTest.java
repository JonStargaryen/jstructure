package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Test;

public class SecondaryStructureElementTest {
    @Test
    public void shouldHandleStructureWithNonStandardAminoAcids() {
        String pdbId = "5tcx";
        Structure structure = StructureParser.source(pdbId).parse();
        structure.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid));
    }

    @Test(expected = ComputationException.class)
    public void shouldFailOnStructureWithNonStandardAminoAcids() {
        String pdbId = "5tcx";
        Structure structure = StructureParser.source(pdbId)
                // skip ligand parsing, thus there is no way to decide whether non-standard groups are actually amino acids
                .minimalParsing(true)
                .parse();
        structure.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid));
    }
}