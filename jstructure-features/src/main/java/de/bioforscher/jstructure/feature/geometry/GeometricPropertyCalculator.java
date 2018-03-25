package de.bioforscher.jstructure.feature.geometry;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;

/**
 * Reports basic geometric properties such as the distance to the centroid and center of mass of individual residues.
 */
public class GeometricPropertyCalculator extends FeatureProvider {
    @Override
    protected void processInternally(Structure structure) {
        structure.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra centerOfMass = chain.calculate().centerOfMass();
        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra centroid = chain.calculate().centroid();

        chain.groups()
                .forEach(group -> {
                    group.getFeatureContainer().addFeature(new GeometricProperties(this,
                            group.calculate().centerOfMass().distance(centerOfMass),
                            group.calculate().centroid().distance(centroid)));
                });
    }
}
