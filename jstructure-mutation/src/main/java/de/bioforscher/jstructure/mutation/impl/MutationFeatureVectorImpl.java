package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.MutationFeatureVector;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.lang.reflect.Method;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

/**
 * The implementation of the data structure to describe the effect of 1 particular amino acid mutation.
 * Created by bittrich on 7/13/17.
 */
public class MutationFeatureVectorImpl implements MutationFeatureVector {
    private static final Logger logger = LoggerFactory.getLogger(MutationFeatureVectorImpl.class);

    private final MutationJob mutationJob;
    private final String mutationIdentifier;
    private final String originalGutteridge;

    private final double exchangeScore;
    private final double evolutionaryInformation;

    private final double maximumAccessibleSurfaceAreaDelta;
    private final double pKaDelta;
    private final boolean gutteridgeChanged;

    private final double functionalGutteridge;

    private final double energy;
    private final double loopFraction;
    private final double rasa;

    private final double energyAminoAcidDelta;
    private final double energyEnvironmentDelta;
    private final double loopFractionAminoAcidDelta;
    private final double loopFractionEnvironmentDelta;
    private final double rasaAminoAcidDelta;
    private final double rasaEnvironmentDelta;

    public MutationFeatureVectorImpl(MutationJob mutationJob,
                                     ChainIdentifier chainIdentifier,
                                     ResidueIdentifier residueIdentifierToMutate,

                                     AminoAcid originalAminoAcid,
                                     AminoAcid mutatedAminoAcid,

                                     List<AminoAcid> originalEnvironment,
                                     List<AminoAcid> mutatedEnvironment,

                                     MutationEffectPredictionServiceImpl.PhysicochemicalFeatureVector fvOriginalAminoAcid,

                                     MutationEffectPredictionServiceImpl.DeltaPhysicochemicalFeatureVector dfvAminoAcid,
                                     MutationEffectPredictionServiceImpl.DeltaPhysicochemicalFeatureVector dfvEnvironment) {
        // general information
        this.mutationJob = mutationJob;
        this.mutationIdentifier = originalAminoAcid.getOneLetterCode() +
                residueIdentifierToMutate.getResidueNumber() +
                mutatedAminoAcid.getOneLetterCode();
        this.originalGutteridge = originalAminoAcid.getGroupPrototype().getGutteridgeGrouping().name();

        // PSSM related values
        Optional<EvolutionaryInformation> featureOptional = originalAminoAcid.getFeatureContainer().getFeatureOptional(EvolutionaryInformation.class);
        this.exchangeScore = featureOptional.map(EvolutionaryInformation::getExchangeScores)
                .map(scores -> scores.get(AminoAcid.Family.resolveGroupPrototype(originalAminoAcid.getGroupPrototype())))
                .orElse(0.0);
        this.evolutionaryInformation = featureOptional.map(EvolutionaryInformation::getInformation).orElse(0.0);

        // generic information derived from the nature of the mutation
        GroupPrototype originalPrototype = originalAminoAcid.getGroupPrototype();
        GroupPrototype mutatedPrototype = mutatedAminoAcid.getGroupPrototype();
        this.maximumAccessibleSurfaceAreaDelta = originalPrototype.getMaximumAccessibleSurfaceArea()
                - mutatedPrototype.getMaximumAccessibleSurfaceArea();
        this.pKaDelta = originalPrototype.getIsoelectricPoint()
                - mutatedPrototype.getIsoelectricPoint();
        this.gutteridgeChanged = originalPrototype.getGutteridgeGrouping() != mutatedPrototype.getGutteridgeGrouping();

        // information derived from the sequential environment
        double environmentSize = originalEnvironment.size();

        // fraction of amino acids with functional Gutteridge grouping
        this.functionalGutteridge = originalEnvironment.stream()
                .map(AminoAcid::getGroupPrototype)
                .map(GroupPrototype::getGutteridgeGrouping)
                .filter(gutteridgeGrouping -> gutteridgeGrouping != GroupPrototype.GutteridgeGrouping.NONE)
                .count() / environmentSize;

        // features of the original amino acid
        this.energy = fvOriginalAminoAcid.getEnergy();
        this.loopFraction = fvOriginalAminoAcid.getLoopFraction();
        this.rasa = fvOriginalAminoAcid.getRasa();

        // information derived from the structural environment
        this.energyAminoAcidDelta = dfvAminoAcid.getEnergy();
        this.energyEnvironmentDelta = dfvEnvironment.getEnergy();
        this.loopFractionAminoAcidDelta = dfvAminoAcid.getLoopFraction();
        this.loopFractionEnvironmentDelta = dfvAminoAcid.getLoopFraction();
        this.rasaAminoAcidDelta = dfvAminoAcid.getRasa();
        this.rasaEnvironmentDelta = dfvEnvironment.getRasa();
    }

    public double[] toDoubleArray() {
        return Stream.of(getClass().getDeclaredMethods())
                .filter(method -> method.getName().startsWith("get") || method.getName().startsWith("is"))
                .filter(method -> !"getClass".equals(method.getName()) &&
                        !"getMutationJob".equals(method.getName()) &&
                        !"getMutationIdentifier".equals(method.getName()))
                .mapToDouble(this::invokeGetterSafely)
                .toArray();
    }

    @Override
    public String toDoubleString() {
        return DoubleStream.of(toDoubleArray())
                .mapToObj(StandardFormat::format)
                .collect(Collectors.joining(","));
    }

    public static String toPartialHeader() {
        return Stream.of(MutationFeatureVectorImpl.class.getDeclaredMethods())
                .filter(method -> method.getName().startsWith("get") || method.getName().startsWith("is"))
                .filter(method -> !"getClass".equals(method.getName()) &&
                        !"getMutationJob".equals(method.getName()) &&
                        !"getMutationIdentifier".equals(method.getName()))
                .map(Method::getName)
                .map(fieldName -> fieldName.replace("get", ""))
                .map(fieldName -> "@ATTRIBUTE " + fieldName + " NUMERIC")
                .collect(Collectors.joining(System.lineSeparator()));
    }

    @Override
    public String toString() {
        return mutationJob.getJobName() + "," +
                mutationIdentifier + "," +
                originalGutteridge + "," +
                toDoubleString();

    }

    private double invokeGetterSafely(Method method) {
        try {
            Object value = method.invoke(this);
            if(value instanceof Double) {
                return (double) value;
            }
            if(value instanceof Boolean) {
                return (Boolean) value ? 1.0 : 0.0;
            }
            throw new IllegalArgumentException("unexpected content type " + value.getClass() + " of method " +
                    method.getName());
        } catch (Exception e) {
            logger.warn("[{}] could not invoke getter method {}",
                    mutationJob.getUuid(),
                    method.getName(),
                    e);
            return 0.0;
        }
    }

    public MutationJob getMutationJob() {
        return mutationJob;
    }

    @Override
    public String getMutationIdentifier() {
        return mutationIdentifier;
    }

    public double getExchangeScore() {
        return exchangeScore;
    }

    public double getEvolutionaryInformation() {
        return evolutionaryInformation;
    }

    public double getMaximumAccessibleSurfaceAreaDelta() {
        return maximumAccessibleSurfaceAreaDelta;
    }

    public boolean isGutteridgeChanged() {
        return gutteridgeChanged;
    }

    public double getFunctionalGutteridge() {
        return functionalGutteridge;
    }

    public double getpKaDelta() {
        return pKaDelta;
    }

    public double getEnergy() {
        return energy;
    }

    public double getLoopFraction() {
        return loopFraction;
    }

    public double getRasa() {
        return rasa;
    }

    public double getEnergyAminoAcidDelta() {
        return energyAminoAcidDelta;
    }

    public double getEnergyEnvironmentDelta() {
        return energyEnvironmentDelta;
    }

    public double getLoopFractionAminoAcidDelta() {
        return loopFractionAminoAcidDelta;
    }

    public double getLoopFractionEnvironmentDelta() {
        return loopFractionEnvironmentDelta;
    }

    public double getRasaAminoAcidDelta() {
        return rasaAminoAcidDelta;
    }

    public double getRasaEnvironmentDelta() {
        return rasaEnvironmentDelta;
    }
}
