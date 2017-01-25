package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.alignment.AbstractAlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Aligns two energy profiles and scores the resulting alignment by the dScore.
 * Created by S on 20.01.2017.
 * @author originally written by Florian Heinke
 */
public class EnergyProfileAligner extends AbstractAlignmentAlgorithm {
    private static final Logger logger = LoggerFactory.getLogger(EnergyProfileAligner.class);
    private static final int DEFAULT_GAP_COST = -40;
    private final int gapCost;
    private static final int DEFAULT_GAP_EXTEND = -7;
    private final int gapExtend;
    private static final int DEFAULT_NUMBER_OF_PERMUTATIONS = 50;
    private final int numberOfPermutations;
    private boolean isSemiGlobal;

    private int[][] pathMatrix;
    private int[][] fMatrix;

    private List<Integer> ei;
    private List<Integer> ej;

    private int scoreReal;

    private static final double[] ENERGY_BOUNDARIES = { -Double.MAX_VALUE, -32.19657709214838, -28.360766398274798, -25.23929221080452, -22.328373156958044, -19.493175687527987, -16.62469913138201, -13.850806967402914,
            -11.369088398962608, -9.314490868789694, -7.6492841662922855, -6.298501258749379, -5.141137053601808, -4.1008253531444865, -3.1690257375033464, -2.27957855620535, -1.389078754253618,
            -0.4313635401829649, 0.7109572780383445, 8.555324534814812, Double.MAX_VALUE };

    public EnergyProfileAligner() {
        this(DEFAULT_GAP_COST, DEFAULT_GAP_EXTEND, DEFAULT_NUMBER_OF_PERMUTATIONS);
    }

    public EnergyProfileAligner(int gapCost, int gapExtend, int numberOfPermutations) {
        this.gapCost = gapCost;
        this.gapExtend = gapExtend;
        this.numberOfPermutations = numberOfPermutations;
        this.isSemiGlobal = false;
    }

    @Override
    public AlignmentResult align(AtomContainer reference, AtomContainer query) {
        List<Integer> referenceEnergyProfile = extractDiscreteEnergyProfile(reference);
        List<Integer> queryEnergyProfile = extractDiscreteEnergyProfile(query);

        if(referenceEnergyProfile.size() > queryEnergyProfile.size()) {
            ei = referenceEnergyProfile;
            ej = queryEnergyProfile;
        } else {
            ei = queryEnergyProfile;
            ej = referenceEnergyProfile;
        }

        //TODO need clone?

        createFMatrix(ei.size(), ej.size());

        scoreReal = calculateAlignment();

        // permutate alignment to obtain significance score
        double distanceScore = calculateDistanceScore();

        return new AlignmentResult(reference, query, null, distanceScore, null, null);
    }

    private double calculateDistanceScore() {
        double sum = 0;
        for (int i = 0; i < numberOfPermutations; i++) {
            List<Integer> permutationList = new ArrayList<>();
            for (int permutation = 0; permutation < ej.size(); permutation++) {
                permutationList.add(permutation);
            }
            int k = 0;
            while (permutationList.size() > 0) {
                int pos = (int) (Math.random() * permutationList.size());
                permutationList.remove(pos);
                int e_1 = ej.get(k);
                int e_2 = ej.get(pos);
                ej.set(pos, e_1);
                ej.set(k, e_2);
                k++;
            }

            permutationList = new ArrayList<>();
            for (int per = 0; per < ei.size(); per++) {
                permutationList.add(per);
            }
            k = 0;
            while (permutationList.size() > 0) {
                int pos = (int) (Math.random() * permutationList.size());
                permutationList.remove(pos);
                int e_1 = ei.get(k);
                int e_2 = ei.get(pos);
                ei.set(pos, e_1);
                ei.set(k, e_2);
                k++;
            }

            sum += calculateAlignment();
        }
        double mean = sum / numberOfPermutations;

        double distanceScore = -Math.log((scoreReal - mean) / (((ei.size() * 12 + ej.size() * 12) / 2) - mean));
        if (Double.isNaN(distanceScore) || distanceScore > 5) {
            distanceScore = 5;
        }
        return distanceScore;
    }

    private int calculateAlignment() {
        for (int i = 1; i < ei.size() + 1; i++) {
            int eiEntry = ei.get(i - 1);
            for (int j = 1; j < ej.size() + 1; j++) {

                int qdiag, qtop, qleft;
                int ejEntry = ej.get(j - 1);

                qdiag = fMatrix[i - 1][j - 1] + calcEScore(eiEntry, ejEntry);

                if (pathMatrix[i - 1][j] != 2) {
                    qtop = fMatrix[i - 1][j] + gapExtend;
                } else {
                    qtop = fMatrix[i - 1][j] + gapCost;
                }

                if (pathMatrix[i][j - 1] != 2) {
                    qleft = fMatrix[i][j - 1] + gapExtend;
                } else {
                    qleft = fMatrix[i][j - 1] + gapCost;
                }

                int tempMax = Math.max(qdiag, qtop);
                int max = Math.max(tempMax, qleft);

                fMatrix[i][j] = max;

                if (max == qdiag) {
                    pathMatrix[i][j] = 2;
                } else if (max == qtop) {
                    pathMatrix[i][j] = 1;
                } else {
                    pathMatrix[i][j] = -1;
                }
            }
        }
        return fMatrix[ei.size()][ej.size()];
    }

    private int calcEScore(int e_i, int e_j) {
        switch (Math.abs(e_i - e_j)) {
            case 0:
                return 12;
            case 1:
                return 8;
            case 2:
                return 6;
            case 3:
                return 2;
            case 4:
                return -4;
            case 5:
                return -8;
            default:
                return -10;
        }
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
            if(ENERGY_BOUNDARIES[index] < solvationEnergy && ENERGY_BOUNDARIES[index + 1] > solvationEnergy) {
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