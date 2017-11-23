package de.bioforscher.equant.train;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.regex.Pattern;

public class TargetFunctionParser extends FeatureProvider {
    public void process(Structure structure, Path resultPath) {
        process(structure, EquantConstants.newInputStream(resultPath));
    }

    private static final Pattern PATTERN = Pattern.compile("\\s+");

    void process(Structure structure, InputStream inputStream) {
        //TODO breaks when selecting directly on structure with NPE - why?
        Chain chain = structure.chains().findFirst().get();
        // molecule 1 is prediction, 2 is target
        // #      Molecule1      Molecule2  DISTANCE    Mis    MC     All    Dist_max   GDC_mc  GDC_all
        // LGA    M       1      M       1    14.220     0    0.535   0.905    15.441    0.000    0.000
        new BufferedReader(new InputStreamReader(inputStream))
                .lines()
                .filter(line -> line.startsWith("LGA "))
                .map(PATTERN::split)
                .forEach(split -> {
                    int residueNumber = Integer.valueOf(split[2]);
                    double distance = Double.valueOf(split[5]);
                    AminoAcid aminoAcid = chain.select()
                            .residueNumber(residueNumber)
                            .asAminoAcid();
                    aminoAcid.getFeatureContainer().addFeature(new TargetFunction(this,
                            distance));
                });
    }

    @Override
    protected void processInternally(Structure structure) {
        throw new UnsupportedOperationException("use process(Structure, Path) instead");
    }

    public static class TargetFunction extends FeatureContainerEntry {
        private static final double DISTANCE_THRESHOLD = 3.0;
        private final double distance;
        private final double sscore;

        public TargetFunction(FeatureProvider featureProvider, double distance) {
            super(featureProvider);
            this.distance = distance;
            double fraction = distance / DISTANCE_THRESHOLD;
            this.sscore = 1 / (1 + fraction * fraction);
        }

        public double getDistance() {
            return distance;
        }

        public double getSscore() {
            return sscore;
        }
    }
}
