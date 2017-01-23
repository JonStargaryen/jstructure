package de.bioforscher.jstructure.feature.loopfraction;

import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Fragment;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Computes the loop fraction of a protein.
 * Created by S on 16.01.2017.
 */
@FeatureProvider(provides = LoopFractionCalculator.LOOP_FRACTION, requires = SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES)
public class LoopFractionCalculator extends AbstractFeatureProvider {
    private static final String RAW_LOOP_FRACTION = "RAW_LOOP_FRACTION";
    public static final String LOOP_FRACTION = "LOOP_FRACTION";
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
    protected void processInternally(Protein protein) {
        assignRawLoopFraction(protein);
        smoothLoopFraction(protein);
    }

    /**
     * Transforms the binary annotation of secondary structure elements into the loop fraction feature by averaging the
     * value with sequential neighbors in a window of size 9.
     * @param protein the container to process
     */
    private void smoothLoopFraction(Protein protein) {
        protein.chains().forEach(this::smoothLoopFraction);
    }

    private void smoothLoopFraction(Chain chain) {
        // smooth values for all groups with the complete set of neighbors
        Combinatorics.fragmentsOf(chain.aminoAcids().collect(Collectors.toList()), WINDOW_SIZE).forEach(this::smoothLoopFraction);

        // smooth values of start and end residues of the chain
        chain.aminoAcids()
                .filter(group -> !group.getFeatureMap().containsKey(LOOP_FRACTION))
                .forEach(group -> {
                    int residueNumber = group.getResidueNumber();
                    double smoothedValue = Selection.on(chain)
                            .residueNumber(surroundingResidueNumbers(residueNumber))
                            .asFilteredGroups()
                            .mapToDouble(surroundingGroup -> surroundingGroup.getFeatureAsDouble(RAW_LOOP_FRACTION))
                            .average()
                            .orElse(1.0);
                    group.setFeature(LOOP_FRACTION, smoothedValue);
        });
    }

    private int[] surroundingResidueNumbers(int residueNumber) {
        int[] surroundingResidueNumbers = new int[WINDOW_SIZE];

        for(int i = 0; i < WINDOW_SIZE; i++) {
            surroundingResidueNumbers[i] = residueNumber - (WINDOW_SIZE / 2) + i;
        }

        return surroundingResidueNumbers;
    }

    private void smoothLoopFraction(Fragment<Group> fragment) {
        Group middleGroup = fragment.getElement(WINDOW_MIDDLE_INDEX);
        double smoothedLoopFraction = 0;

        for(int i = 0; i < WINDOW_SIZE; i++) {
            smoothedLoopFraction += weights[i] * fragment.getElement(i).getFeatureAsDouble(RAW_LOOP_FRACTION);
        }
        middleGroup.setFeature(LOOP_FRACTION, smoothedLoopFraction);
    }

    /**
     * Assigns the raw loop fraction, i.e. 1 if loop and 0 otherwise.
     * @param protein the container to process
     */
    private void assignRawLoopFraction(Protein protein) {
        protein.aminoAcids().forEach(this::assignRawLoopFraction);
    }

    private void assignRawLoopFraction(Group group) {
        group.setFeature(RAW_LOOP_FRACTION, group.getFeature(SecStrucState.class,
                    SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).getSecondaryStructure().isCoilType() ? 1.0 : 0.0);
    }
}
