package de.bioforscher.jstructure.feature.cerosene;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.stream.IntStream;

/**
 * An adaptation of the sequence cerosene tool which transforms 3D coordinates of a protein in the RGB color space. This
 * allows the identification of residues in spatial proximity on sequence level as they will feature identical color
 * values.
 * Created by bittrich on 3/14/17.
 * originally proposed by Florian Heinke
 * SequenceCEROSENE: a computational method and web server to visualize spatial residue neighborhoods at the sequence level.
 * Florian Heinke, Sebastian Bittrich, Florian Kaiser and Dirk Labudde. BioData Mining 2016 9:6. 10.1186/s13040-016-0083-7
 */
@FeatureProvider(provides = { CeroseneEmbedding.class })
public class SequenceCerosene extends AbstractFeatureProvider {
    private static final double[] SUMMAND = { 127, 127, 127 };

    @Override
    protected void processInternally(Protein protein) {
        final double[] centroid = LinearAlgebraAtom.centroid(protein);
        final double maximalExtent = LinearAlgebra3D.capToInterval(23.0, maximalExtent(protein, centroid), Double.MAX_VALUE);

        protein.chains()
                .flatMap(Chain::groups)
                .forEach(group -> {
                    double[] coordinates = LinearAlgebraAtom.centroid(group);
                    double[] rgb = LinearAlgebra3D.add(LinearAlgebra3D.multiply(LinearAlgebra3D.subtract(centroid, coordinates), (double) 127 / maximalExtent), SUMMAND);
                    double[] hsv = rgb2hsv(rgb);
                    group.getFeatureContainer().addFeature(new CeroseneEmbedding(this, rgb, hsv));
                });
    }

    static double maximalExtent(Protein protein, double[] centroid) {
        return protein.atoms()
                .map(Atom::getCoordinates)
                .map(coordinates -> LinearAlgebra3D.subtract(coordinates, centroid))
                .mapToDouble(difference -> IntStream.range(0, 3).mapToDouble(dimension -> difference[dimension]).max().getAsDouble())
                .max()
                .orElse(23.0);
    }


    static double[] rgb2hsv(double[] rgb) {
        rgb = LinearAlgebra3D.divide(rgb, 255.0);

        double min = Math.min(Math.min(rgb[0], rgb[1]), rgb[2]);
        double max = Math.max(Math.max(rgb[0], rgb[1]), rgb[2]);
        double delta = max - min;

        double h = max;
        double s;

        if(delta == 0) {
            h = 0;
            s = 0;
        } else {
            s = delta / max;

            double delR = ( ( ( max - rgb[0] ) / 6 ) + ( delta / 2 ) ) / delta;
            double delG = ( ( ( max - rgb[1] ) / 6 ) + ( delta / 2 ) ) / delta;
            double delB = ( ( ( max - rgb[2] ) / 6 ) + ( delta / 2 ) ) / delta;

            if(rgb[0] == max) {
                h = delB - delG;
            } else if(rgb[1] == max) {
                h = (1.0 / 3.0) + delR - delB;
            } else if(rgb[2] == max) {
                h = (2.0 / 3.0) + delG - delR;
            }

            if(h < 0) {
                h += 1;
            }
            if(h > 1) {
                h -= 1;
            }
        }

        return new double[] { h, s, max };
    }
}
