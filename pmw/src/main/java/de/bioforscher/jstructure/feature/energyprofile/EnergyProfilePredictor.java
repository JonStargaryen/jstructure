package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * An adapation of the GOR-algorithm to predict energy profiles from a protein sequence.
 * Created by S on 20.01.2017.
 * @author originally written by Florian Heinke
 */
@FeatureProvider(provides = EnergyProfileCalculator.SOLVATION_ENERGY)
public class EnergyProfilePredictor extends AbstractFeatureProvider {
    private final int windowSize;
    private final String spacer;

    private static final int DEFAULT_WINDOW_SIZE = 3;
    private static final double[] QUANTIL_ENERGIES = { -1.17, -5.24, -13.3, -27.33 };

    public EnergyProfilePredictor() {
        this(DEFAULT_WINDOW_SIZE);
    }

    public EnergyProfilePredictor(int windowSize) {
        this.windowSize = windowSize;
        this.spacer = IntStream.range(0, windowSize).mapToObj(i -> "#").collect(Collectors.joining());
        initializeLibrary();
    }

    private void initializeLibrary() {

    }

    @Override
    protected void processInternally(Protein protein) {
        if(protein.getAminoAcidSequence().length() < 7) {
            throw new IllegalArgumentException("cannot predict energy profile for sequences shorter than 7 residues");
        }

        String sequence = spacer + protein.getAminoAcidSequence() + spacer;
    }
}
