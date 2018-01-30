package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.mathematics.Fragment;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Computes the loop fraction of a protein.
 * Created by S on 16.01.2017.
 */
public class LoopFractionCalculator extends FeatureProvider {
    private static final int WINDOW_SIZE = 9;
    private static final int WINDOW_MIDDLE_INDEX = (WINDOW_SIZE / 2) + 1;
    private final double[] weights;

    public LoopFractionCalculator() {
        weights = weights(WINDOW_SIZE);
    }

    /**
     * Computes the gaussian weights for a bell-shaped smoothing rule.
     * @param windowSize the window size
     * @return the computed weights
     */
    private double[] weights(int windowSize) {
        int halfWindowSize = windowSize / 2;
        double sigma = Math.sqrt(halfWindowSize);
        double[] weights = new double[windowSize];
        for(int i = 0; i < windowSize; i++) {
            weights[i] = 1 / (sigma * Math.sqrt(2 * Math.PI)) * Math.pow(Math.E, -0.5 * Math.pow((i - halfWindowSize) / sigma, 2));
        }

        double sum = 1 - DoubleStream.of(weights).sum();
        for(int i = 0; i < windowSize; i++) {
            weights[i] += sum / windowSize;
        }

        return weights;
    }

    @Override
    protected void processInternally(Structure protein) {
        assignRawLoopFraction(protein);
        smoothLoopFraction(protein);
    }

    /**
     * Transforms the binary annotation of secondary structure elements into the loop fraction feature by averaging the
     * value with sequential neighbors in a window of size 9.
     * @param protein the container to processUniProtId
     */
    private void smoothLoopFraction(Structure protein) {
        protein.chains()
                .filter(chain -> chain.aminoAcids().count() > WINDOW_SIZE)
                .forEach(this::smoothLoopFraction);
    }

    private void smoothLoopFraction(Chain chain) {
        // smooth values for all groups with the complete set of neighbors
        SetOperations.fragmentsOf(chain.aminoAcids()
                        .collect(Collectors.toList()), WINDOW_SIZE)
                .forEach(this::smoothLoopFraction);

        // smooth values of start and end residues of the chain
        chain.aminoAcids()
                .filter(group -> !group.getFeatureContainer().getFeatureOptional(LoopFraction.class).isPresent())
                .forEach(group -> {
                    int residueNumber = group.getResidueIdentifier().getResidueNumber();
                    double smoothedValue = chain.select()
                            .aminoAcids()
                            .residueNumber(surroundingResidueNumbers(residueNumber))
                            .asFilteredGroups()
                            .mapToDouble(surroundingGroup -> surroundingGroup.getFeature(RawLoopFraction.class).getRawLoopFraction())
                            .average()
                            .orElse(1.0);
                    group.getFeatureContainer().addFeature(new LoopFraction(this, smoothedValue));
        });
    }

    private int[] surroundingResidueNumbers(int residueNumber) {
        int[] surroundingResidueNumbers = new int[WINDOW_SIZE];

        for(int i = 0; i < WINDOW_SIZE; i++) {
            surroundingResidueNumbers[i] = residueNumber - (WINDOW_SIZE / 2) + i;
        }

        return surroundingResidueNumbers;
    }

    private void smoothLoopFraction(Fragment<AminoAcid> fragment) {
        Group middleGroup = fragment.getElement(WINDOW_MIDDLE_INDEX);
        double smoothedLoopFraction = 0;

        for(int i = 0; i < WINDOW_SIZE; i++) {
            smoothedLoopFraction += weights[i] * fragment.getElement(i).getFeature(RawLoopFraction.class).getRawLoopFraction();
        }
        middleGroup.getFeatureContainer().addFeature(new LoopFraction(this, smoothedLoopFraction));
    }

    /**
     * Assigns the raw loop fraction, i.e. 1 if loop and 0 otherwise.
     * @param protein the container to processUniProtId
     */
    private void assignRawLoopFraction(Structure protein) {
        protein.aminoAcids()
                .forEach(this::assignRawLoopFraction);
    }

    private void assignRawLoopFraction(Group group) {
        group.getFeatureContainer().addFeature(new RawLoopFraction(this,
                group.getFeature(DSSPSecondaryStructure.class).getSecondaryStructure().isCoilType() ? 1.0 : 0.0));
    }
}
