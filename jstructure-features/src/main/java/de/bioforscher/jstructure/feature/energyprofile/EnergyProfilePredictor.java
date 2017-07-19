package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * An adaptation of the GOR-algorithm to predict energy profiles from a protein sequence.
 * Created by S on 20.01.2017.
 * @author originally written by Florian Heinke
 */
@FeatureProvider(provides = EnergyProfile.class, origin = FeatureProvider.FeatureOrigin.PREDICTION, priority = 200)
public class EnergyProfilePredictor extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(EnergyProfilePredictor.class);
    private static final int WINDOW_SIZE = 3;
    public static final double[] QUANTIL_BOUNDARIES = { -3.17, -8.10, -19.65 };
    public static final double[] QUANTIL_ENERGIES = { -1.17, -5.24, -13.3, -27.33 };
    private static final String BASE_PATH = "energyprofile/";
    private static final String OCCURRENCE_GLOBULAR_PATH = BASE_PATH + "egor-occurrence-globular.dat";
    private static final String PAIRINGS_GLOBULAR_PATH = BASE_PATH + "egor-pairings-globular.dat";
    private Map<String, Map<Integer, Integer>> OCCURRENCE_GLOBULAR;
    // example entry: AA_-3|1 -> 47
    private Map<String, Integer> PAIRINGS_GLOBULAR;
    private final String spacer;

    public EnergyProfilePredictor() {
        initializeLibrary();

        this.spacer = IntStream.range(0, WINDOW_SIZE)
                .mapToObj(i -> "#")
                .collect(Collectors.joining());
    }

    private void initializeLibrary() {
        // load occurrence data
        OCCURRENCE_GLOBULAR = getResourceAsStream(OCCURRENCE_GLOBULAR_PATH)
                .filter(line -> !line.startsWith("\t"))
                .map(line -> line.split("\t"))
                .collect(Collectors.toMap(split -> split[0], this::composeOccurrenceMap));

        // load pairings data
        PAIRINGS_GLOBULAR = getResourceAsStream(PAIRINGS_GLOBULAR_PATH)
                .map(line -> line.split("\t"))
                .collect(Collectors.toMap(line -> line[0], line -> Integer.valueOf(line[1])));
    }

    private Map<Integer, Integer> composeOccurrenceMap(String[] split) {
        return IntStream.range(1, 6)
                .boxed()
                .collect(Collectors.toMap(Function.identity(), i -> Integer.valueOf(split[i])));
    }

    public static double mapEnergyToQuantil(double energy) {
        if (energy > QUANTIL_BOUNDARIES[0]) {
            return QUANTIL_ENERGIES[0];
        } else if(energy > QUANTIL_BOUNDARIES[1]) {
            return QUANTIL_ENERGIES[1];
        } else if(energy > QUANTIL_BOUNDARIES[2]) {
            return QUANTIL_ENERGIES[2];
        } else {
            return QUANTIL_ENERGIES[3];
        }
    }

    @Override
    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids()
                .forEach(this::predictEnergyProfile);
    }

    private void predictEnergyProfile(Chain chain) {
        if(chain.getAminoAcidSequence().length() < 7) {
            logger.warn("cannot predict energy profile for chains with less than 7 residues - found sequence {}", chain.getAminoAcidSequence());
            return;
        }

        String sequence = spacer + chain.getAminoAcidSequence() + spacer;
        List<QuantileEntropySum[]> entropies = calculateEntropies(sequence);
        predictByEntropy(chain, entropies);
    }

    public List<Integer> predictEnergyProfile(String sequence) {
        return calculateEntropies(spacer + sequence + spacer).stream()
                .map(this::getMaximalEntropy)
                .map(maximalEntropy -> maximalEntropy.getIndex() - 1)
                .collect(Collectors.toList());
    }

    private void predictByEntropy(Chain chain, List<QuantileEntropySum[]> entropies) {
        List<Group> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        for(int i = 0; i < entropies.size(); i++) {
            QuantileEntropySum maximalEntropy = getMaximalEntropy(entropies.get(i));
            aminoAcids.get(i)
                    .getFeatureContainer()
                    .addFeature(new EnergyProfile(this, QUANTIL_ENERGIES[maximalEntropy.getIndex() - 1]));
        }
    }

    private QuantileEntropySum getMaximalEntropy(QuantileEntropySum[] entropySums) {
        return Stream.of(entropySums)
                .max(Comparator.comparingDouble(QuantileEntropySum::getEntropy))
                .orElseThrow(() -> new IllegalArgumentException("found no maximal entropy value"));
    }

    private List<QuantileEntropySum[]> calculateEntropies(String sequenceWithSpacers) {
        List<QuantileEntropySum[]> sums = new ArrayList<>();

        for(int i = WINDOW_SIZE; i < sequenceWithSpacers.length() - WINDOW_SIZE; i++) {
            String aa_i = sequenceWithSpacers.substring(i, i + 1);
            QuantileEntropySum[] results = new QuantileEntropySum[4];

            for(int quantile = 1; quantile < 5; quantile++) {
                int c = -WINDOW_SIZE;

                double occurrence = OCCURRENCE_GLOBULAR.get(aa_i).get(quantile);
                double n_occurrence = OCCURRENCE_GLOBULAR.get(aa_i).get(5) - occurrence;
                double addend = Math.log(occurrence / n_occurrence);

                double sum = 0;
                for(int j = i - WINDOW_SIZE; j < i + WINDOW_SIZE + 1; j++) {
                    if(c == 0) {
                        continue;
                    }

                    String aa_j = sequenceWithSpacers.substring(j, j + 1);
                    // an entry is e.g. AA_-3|1
                    double joint_occurrence = PAIRINGS_GLOBULAR.getOrDefault(aa_i + aa_j + "_" + c + "|" + quantile, 1);

                    double n_joint_occurrence = 0.0;
                    for(int counter_quantile = 1; counter_quantile < 5; counter_quantile++) {
                        if(counter_quantile == quantile) {
                            continue;
                        }

                        n_joint_occurrence += PAIRINGS_GLOBULAR.getOrDefault(aa_i + aa_j + "_" + c + "|" + counter_quantile, 1);
                    }

                    sum += Math.log(joint_occurrence / n_joint_occurrence) + Math.log(n_occurrence / occurrence);
                    c++;
                }

                double entropy = addend + sum;
                results[quantile - 1] = new QuantileEntropySum(entropy, quantile, aa_i);
            }

            sums.add(results);
        }

        return sums;
    }

    static class QuantileEntropySum {
        double entropy;
        int index;
        String aminoAcid;

        QuantileEntropySum(double entropy, int index, String aminoAcid) {
            this.entropy = entropy;
            this.index = index;
            this.aminoAcid = aminoAcid;
        }

        double getEntropy() {
            return entropy;
        }

        int getIndex() {
            return index;
        }
    }
}
