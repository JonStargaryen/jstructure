package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.feature.sse.SecondaryStructureType;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * The functional test for the loop fraction calculator.
 * Created by bittrich on 1/18/17.
 */
public class LoopFractionCalculatorTest {
    private Structure protein;
    private FeatureProvider featureProvider;

    @Before
    public void setup() {
        protein = StructureParser.source("1acj").parse();
        featureProvider = new LoopFractionCalculator();
    }

    @Test
    public void shouldComputeLoopFraction() {
        featureProvider.process(protein);
        protein.aminoAcids().forEach(group -> {
            SecondaryStructureType state = group.getFeature(DSSPSecondaryStructure.class).getSecondaryStructure();
            boolean isCoil = state.isCoilType();
            double value = group.getFeature(LoopFraction.class).getLoopFraction();
            System.out.println(group.getIdentifier() + "\t" + state.getOneLetterRepresentation() + "\t" + value);
            if(!isCoil) {
                Assert.assertTrue(value < 1.0);
            }
        });
    }
}