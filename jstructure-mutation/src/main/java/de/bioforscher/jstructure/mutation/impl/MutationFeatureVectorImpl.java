package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtFeatureContainer;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.ConservationProfile;
import de.bioforscher.jstructure.mutation.MutationFeatureVector;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

import java.lang.reflect.Method;
import java.util.List;
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

    private final double sequenceConservationScore;
    private final double energyConservationScore;

    private final double maximumAccessibleSurfaceAreaDelta;
    private final double pKaDelta;
    private final boolean gutteridgeChanged;

    private final double functionalGutteridge;
    private final double bindingSiteFrequency;
    private final double mutagenFrequency;
    private final double variantFrequency;
    private final double disulfidFrequency;

    private final double energyAminoAcidDelta;
    private final double energyEnvironmentDelta;
    private final double ligandContactsAminoAcidDelta;
    private final double ligandContactsEnvironmentDelta;
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

                                     MutationEffectPredictionServiceImpl.DeltaPhysicochemicalFeatureVector dfvAminoAcid,
                                     MutationEffectPredictionServiceImpl.DeltaPhysicochemicalFeatureVector dfvEnvironment) {
        // general information
        this.mutationJob = mutationJob;
        this.mutationIdentifier = chainIdentifier.getChainId() + "_" +
                originalAminoAcid.getOneLetterCode() +
                residueIdentifierToMutate.getResidueNumber() +
                mutatedAminoAcid.getOneLetterCode();

        // global information such as conservation scores
        ConservationProfile conservationProfile = originalAminoAcid.getFeature(ConservationProfile.class);
        this.sequenceConservationScore = conservationProfile.getSequentialConservation();
        this.energyConservationScore = conservationProfile.getEnergeticalConservation();

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
        List<UniProtFeatureContainer> uniProtFeatureContainers = originalEnvironment.stream()
                .map(aminoAcid -> aminoAcid.getFeature(UniProtFeatureContainer.class))
                .collect(Collectors.toList());

        // fraction of amino acids with functional Gutteridge grouping
        this.functionalGutteridge = originalEnvironment.stream()
                .map(AminoAcid::getGroupPrototype)
                .map(GroupPrototype::getGutteridgeGrouping)
                .filter(gutteridgeGrouping -> gutteridgeGrouping != GroupPrototype.GutteridgeGrouping.NONE)
                .count() / environmentSize;

        // fraction of how often certain UniProt features where annotated for homologous - normalized by #groups * #homologousSequences
        double uniProtNormalization = environmentSize * mutationJob.getHomologousSequences().size();
        this.bindingSiteFrequency = uniProtFeatureContainers.stream()
                .mapToInt(container -> container.getFeatures(FeatureType.CA_BIND,
                        FeatureType.ZN_FING,
                        FeatureType.DNA_BIND,
                        FeatureType.NP_BIND,
                        FeatureType.METAL,
                        FeatureType.BINDING,
                        FeatureType.LIPID).size())
                .sum() / uniProtNormalization;
        this.mutagenFrequency = uniProtFeatureContainers.stream()
                .mapToInt(container -> container.getFeatures(FeatureType.MUTAGEN).size())
                .sum() / uniProtNormalization;
        this.variantFrequency = uniProtFeatureContainers.stream()
                .mapToInt(container -> container.getFeatures(FeatureType.VARIANT).size())
                .sum() / uniProtNormalization;
        this.disulfidFrequency = uniProtFeatureContainers.stream()
                .mapToInt(container -> container.getFeatures(FeatureType.DISULFID).size())
                .sum() / uniProtNormalization;

        // information derived from the structural environment
        this.energyAminoAcidDelta = dfvAminoAcid.getEnergy();
        this.energyEnvironmentDelta = dfvEnvironment.getEnergy();
        this.ligandContactsAminoAcidDelta = dfvAminoAcid.getLigandContacts();
        this.ligandContactsEnvironmentDelta = dfvEnvironment.getLigandContacts();
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
                DoubleStream.of(toDoubleArray()).mapToObj(StandardFormat::format).collect(Collectors.joining(","));

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

    public String getMutationIdentifier() {
        return mutationIdentifier;
    }

    public double getSequenceConservationScore() {
        return sequenceConservationScore;
    }

    public double getEnergyConservationScore() {
        return energyConservationScore;
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

    public double getBindingSiteFrequency() {
        return bindingSiteFrequency;
    }

    public double getMutagenFrequency() {
        return mutagenFrequency;
    }

    public double getVariantFrequency() {
        return variantFrequency;
    }

    public double getDisulfidFrequency() {
        return disulfidFrequency;
    }

    public double getEnergyAminoAcidDelta() {
        return energyAminoAcidDelta;
    }

    public double getEnergyEnvironmentDelta() {
        return energyEnvironmentDelta;
    }

    public double getLigandContactsAminoAcidDelta() {
        return ligandContactsAminoAcidDelta;
    }

    public double getLigandContactsEnvironmentDelta() {
        return ligandContactsEnvironmentDelta;
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
