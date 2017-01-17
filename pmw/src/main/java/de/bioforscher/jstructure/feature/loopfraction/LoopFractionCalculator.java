package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;

/**
 * Computes the loop fraction of a protein.
 * Created by S on 16.01.2017.
 */
@FeatureProvider(providedFeatures = LoopFractionCalculator.LOOP_FRACTION)
public class LoopFractionCalculator extends AbstractFeatureProvider {
    public static final String LOOP_FRACTION = "LOOP_FRACTION";

    @Override
    protected void processInternally(Protein protein) {
        //TODO impl
    }
}
