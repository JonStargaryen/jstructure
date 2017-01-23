package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Aligns two energy profiles and scores the resulting alignment by the dScore.
 * Created by S on 20.01.2017.
 * @author originally written by Florian Heinke
 */
public class EnergyProfileAligner extends AbstractAlignmentAlgorithm {
    private static final Logger logger = LoggerFactory.getLogger(EnergyProfileAligner.class);
    private int scoreReal = -1;
    private int gapCost = -40;
    private int gapExtend = -7;
    private int numberOfPermutations = 50;

    private double distanceScore = -1;
    private double zScore = -1;

    private boolean isSemiGlobal = false;
    private boolean dataGotSwapped = false;
    private boolean isPermutated = false;

    private int[][] pathMatrix;
    private int[][] pathMatrixCopy;
    private int[][] fMatrix;

    private List<Integer> ei;
    private List<Integer> ej;
    private List<Integer> eiCopy;
    private List<Integer> ejCopy;

    private static final double[] ENERGY_BOUNDARIES = { -Double.MAX_VALUE, -32.19657709214838, -28.360766398274798, -25.23929221080452, -22.328373156958044, -19.493175687527987, -16.62469913138201, -13.850806967402914,
            -11.369088398962608, -9.314490868789694, -7.6492841662922855, -6.298501258749379, -5.141137053601808, -4.1008253531444865, -3.1690257375033464, -2.27957855620535, -1.389078754253618,
            -0.4313635401829649, 0.7109572780383445, 8.555324534814812, Double.MAX_VALUE };

    @Override
    public AlignmentResult align(AtomContainer atomContainer1, AtomContainer atomContainer2) {
        List<Integer> energyProfile1 = extractDiscreteEnergyProfile(atomContainer1);
        List<Integer> energyProfile2 = extractDiscreteEnergyProfile(atomContainer2);

        if(energyProfile1.size() > energyProfile2.size()) {
            ei = energyProfile1;
            ej = energyProfile2;
        } else {
            ei = energyProfile2;
            ej = energyProfile1;
            dataGotSwapped = true;
        }

        //TODO need clone?

        createFMatrix(ei.size(), ej.size());

        return null;
    }

    private void createFMatrix(int i, int j) {
        fMatrix = new int[i + 1][j + 1];
        pathMatrix = new int[i + 1][j + 1];

        fMatrix[0][0] = 0;
        if (isSemiGlobal) {
            fMatrix[0][1] = 0;
            fMatrix[1][0] = 0;
        } else {
            fMatrix[0][1] = gapCost;
            fMatrix[1][0] = gapCost;
        }

        pathMatrix[0][0] = 0;

        if (isSemiGlobal) {
            for (int n = 2; n <= j; n++) {
                fMatrix[0][n] = 0;
                pathMatrix[0][n] = -1;
            }
            for (int n = 2; n <= i; n++) {
                fMatrix[n][0] = 0;
                pathMatrix[n][0] = 1;
            }
        } else {
            for (int n = 2; n <= j; n++) {
                fMatrix[0][n] = fMatrix[0][n - 1] + gapExtend;
                pathMatrix[0][n] = -1;
            }
            for (int n = 2; n <= i; n++) {
                fMatrix[n][0] = fMatrix[n - 1][0] + gapExtend;
                pathMatrix[n][0] = 1;
            }
        }
    }

    /**
     * Extracts the discrete energy profile for a given collection of residues.
     * @param container the container to process
     * @return a list of all discretized energy values
     */
    private List<Integer> extractDiscreteEnergyProfile(AtomContainer container) {
        return container.atoms()
                // map to group level - will still collect 'dangling' atoms into a group
                .map(Atom::getParentGroup)
                .distinct()
                // ensure we are only dealing with amino acids
                .filter(Group::isAminoAcid)
                .mapToDouble(this::getSolvationEnergy)
                .mapToInt(this::discretize)
                //TODO is this really faster than working on POJO-stream all the way?
                .boxed()
                .collect(Collectors.toList());
    }

    /**
     * Translates a solvation energy value to its discrete representation.
     * @param solvationEnergy the value
     * @return the interval in which this value is located
     */
    private int discretize(double solvationEnergy) {
        for(int index = 0; index < ENERGY_BOUNDARIES.length - 1; index++) {
            if(ENERGY_BOUNDARIES[index] > solvationEnergy && ENERGY_BOUNDARIES[index + 1] < solvationEnergy) {
                return index;
            }
        }

        throw new IllegalArgumentException("solvation energy value is out of range: " + solvationEnergy);
    }

    /**
     * Merely retrieves the solvation energy value for a given group.
     * @param group the container to process
     * @return its solvation energy value as double
     */
    private double getSolvationEnergy(Group group) {
        return group.getFeatureAsDouble(EnergyProfileCalculator.SOLVATION_ENERGY);
    }
}
