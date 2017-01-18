package feature.loopfraction;

import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.sse.DSSPSecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
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
        protein = ProteinParser.parseProteinById("1acj");
        featureProvider = FeatureProviderRegistry.getInstance().resolve(LoopFractionCalculator.LOOP_FRACTION);
    }

    @Test
    public void shouldComputeLoopFraction() {
        featureProvider.process(protein);
        protein.aminoAcids().forEach(group -> {
            DSSPSecondaryStructureElement state = group.getFeature(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure();
            boolean isCoil = state.isCoilType();
            double value = group.getFeatureAsDouble(LoopFractionCalculator.LOOP_FRACTION);
            System.out.println(group.getIdentifier() + "\t" + state.getOneLetterRepresentation() + "\t" + value);
            if(!isCoil) {
                Assert.assertTrue(value < 1.0);
            }
        });
    }
}