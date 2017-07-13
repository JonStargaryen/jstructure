package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.model.structure.Protein;

/**
 * Annotate needed features in a standardized manner.
 * Created by bittrich on 7/12/17.
 */
class CommonFeatureAnnotator {
    private static final AccessibleSurfaceAreaCalculator accessibleSurfaceAreaCalculator = new AccessibleSurfaceAreaCalculator();
    private static final EnergyProfileCalculator energyProfileCalculator = new EnergyProfileCalculator();
    private static final LoopFractionCalculator loopFractionCalculator = new LoopFractionCalculator();
    private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();

    static void annotateProtein(Protein protein) {
        accessibleSurfaceAreaCalculator.process(protein);
        energyProfileCalculator.process(protein);
        loopFractionCalculator.process(protein);
//        plipAnnotator.process(protein);
    }
}
