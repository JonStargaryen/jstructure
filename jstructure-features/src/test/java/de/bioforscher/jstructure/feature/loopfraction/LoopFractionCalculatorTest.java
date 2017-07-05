package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * The functional test for the loop fraction calculator.
 * Created by bittrich on 1/18/17.
 */
public class LoopFractionCalculatorTest {
    private Protein protein;
    private AbstractFeatureProvider featureProvider;

    @Before
    public void setup() {
        protein = ProteinParser.source("1acj").parse();
        featureProvider = FeatureProviderRegistry.resolve(LoopFraction.class);
    }

    @Test
    public void shouldComputeLoopFraction() {
        featureProvider.process(protein);
        protein.aminoAcids().forEach(group -> {
            SecondaryStructureElement state = group.getFeatureContainer().getFeature(DSSPSecondaryStructure.class).getSecondaryStructure();
            boolean isCoil = state.isCoilType();
            double value = group.getFeatureContainer().getFeature(LoopFraction.class).getLoopFraction();
            System.out.println(group.getIdentifier() + "\t" + state.getOneLetterRepresentation() + "\t" + value);
            if(!isCoil) {
                Assert.assertTrue(value < 1.0);
            }
        });
    }
}