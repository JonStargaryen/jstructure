package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.ConservationCalculator;
import de.bioforscher.jstructure.mutation.ConservationProfile;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Implementation of the energy profile conservation calculator.
 * Created by bittrich on 7/17/17.
 */
public class ConservationCalculatorImpl implements ConservationCalculator {
    private static final Logger logger = LoggerFactory.getLogger(ConservationCalculatorImpl.class);
    private final Rate4SiteWrapper rate4siteWrapper;
    private final EnergyProfileCalculator energyProfileCalculator;

    public ConservationCalculatorImpl() {
        this.rate4siteWrapper = new Rate4SiteWrapper();
        this.energyProfileCalculator = new EnergyProfileCalculator();
    }

    @Override
    public void extractConservationProfile(MutationJob mutationJob) {
        try {
            // write temporary structures
            Path structurePath = Files.createTempDirectory("mmm-in");
            Iterator<Chain> iterator = Stream.concat(Stream.of(mutationJob.getReferenceChain()),
                    mutationJob.getHomologousPdbChains().stream()).iterator();
            while(iterator.hasNext()) {
                Chain chain = iterator.next();
                Files.write(structurePath.resolve(chain.getChainIdentifier().getFullName() + ".pdb"),
                        (chain.getPdbRepresentation()).getBytes());
            }

            // profile is assigned to reference chain
            String sequenceAlignmentString = composeSequenceAlignmentString(mutationJob);
            logger.info("[{}] computing sequence conservation profile by {}",
                    mutationJob.getUuid(),
                    rate4siteWrapper.getClass().getSimpleName());
            List<Double> sequenceScores = rate4siteWrapper.executeCommand(sequenceAlignmentString);

            String energyAlignmentString = composeEnergyAlignmentString(mutationJob);
            logger.info("[{}] computing energy conservation profile by {}",
                    mutationJob.getUuid(),
                    rate4siteWrapper.getClass().getSimpleName());
            List<Double> energyScores = rate4siteWrapper.executeCommand(energyAlignmentString);

            List<AminoAcid> aminoAcids = mutationJob.getReferenceChain().aminoAcids().collect(Collectors.toList());
            if(sequenceScores.size() != aminoAcids.size() || energyScores.size() != aminoAcids.size()) {
                throw new IllegalArgumentException("number of result lines and amino acid count does not match - potentially the alignmentMap is not sorted - the reference sequence should be first");
            }

            for(int i = 0; i < aminoAcids.size(); i++) {
                double sequenceScore = sequenceScores.get(i);
                double energyScore = energyScores.get(i);

                aminoAcids.get(i).getFeatureContainer().addFeature(new ConservationProfile(sequenceScore,
                        energyScore));
            }
        } catch (IOException | InterruptedException e) {
            logger.warn("[{}] encountered exception during conservation profile computation",
                    mutationJob.getUuid(),
                    e);
            throw new ComputationException(e);
        }
    }

    private String composeSequenceAlignmentString(MutationJob mutationJob) {
        return mutationJob.getAlignmentMap().entrySet().stream()
                .filter(entry -> !entry.getKey().contains("_") || entry.getKey().equals(mutationJob.getReferenceChain().getChainIdentifier().getFullName()))
                .map(entry -> ">" + entry.getKey() + System.lineSeparator() + entry.getValue())
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private String composeEnergyAlignmentString(MutationJob mutationJob) {
        return mutationJob.getAlignmentMap().entrySet().stream()
                .filter(entry -> entry.getKey().contains("_"))
                .map(entry -> predictEnergyProfile(mutationJob, entry))
                .map(pair-> ">" + pair.getLeft() + System.lineSeparator() + pair.getRight())
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private Pair<String, String> predictEnergyProfile(MutationJob mutationJob, Map.Entry<String, String> entry) {
        String alignment = entry.getValue();
        Structure protein = Stream.concat(Stream.of(mutationJob.getReferenceChain()), mutationJob.getHomologousPdbChains().stream())
                .filter(chain -> chain.getChainIdentifier().getFullName().equalsIgnoreCase(entry.getKey()))
                .findFirst()
                .get()
                .getParentStructure();
        energyProfileCalculator.process(protein);
        List<String> profile = protein.aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(EnergyProfile.class))
                .map(EnergyProfile::getSolvationEnergy)
                .map(this::mapToDiscreteAlphabet)
                .collect(Collectors.toList());

        // create string representation of discrete energy profile with respect to alignment gaps
        StringBuilder discreteEnergyProfile = new StringBuilder();
        int consumed = 0;
        for(int i = 0; i < alignment.length(); i++) {
            char alignmentChar = alignment.charAt(i);
            if(alignmentChar == '-') {
                discreteEnergyProfile.append(alignmentChar);
            } else {
                discreteEnergyProfile.append(profile.get(consumed));
                consumed++;
            }
        }

        return new Pair<>(entry.getKey(), discreteEnergyProfile.toString());
    }

    /**
     * 4-letter representation of energy quantiles. 'A' is highest energy, 'T' lowest.
     */
    private static final String[] REPRESENTATION = { "A", "C", "G", "T" };

    private String mapToDiscreteAlphabet(int i) {
        return REPRESENTATION[i];
    }

    private String mapToDiscreteAlphabet(double energy) {
//        QUANTIL_ENERGIES = { -1.17, -5.24, -13.3, -27.33 };
        if(energy > 0.5 * (-1.17 + -5.24)) {
            return "A";
        } else if(energy > 0.5 * (-5.24 + -13.3)) {
            return "C";
        } else if(energy > 0.5 * (-13.3 + -27.33)) {
            return "G";
        } else {
            return "T";
        }
    }
}
