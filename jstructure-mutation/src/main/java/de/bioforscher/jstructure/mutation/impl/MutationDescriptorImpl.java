package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.StandardNumberFormat;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtFeatureContainer;
import de.bioforscher.jstructure.model.feature.FeatureContainer;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.MutationDescriptor;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

/**
 * Wraps all gathered information for a mutation.
 * Created by bittrich on 7/12/17.
 */
public class MutationDescriptorImpl implements MutationDescriptor {
    /*
     general information
     */
    private final String sequenceIdentifier;
    private final String mutationDescription;

    /*
     gutteridge-specific data
     */
    private final GroupPrototype.GutteridgeGrouping originalGutteridge;
    private final GroupPrototype.GutteridgeGrouping mutatedGutteridge;
    private final boolean originalAminoAcidFunctional;
    private final boolean mutatedAminoAcidFunctional;
    private final boolean mutationChangedGutteridge;

    /*
     change in physicochemical properties due to the mutation TODO BLOSUM?
     */
    private final double maximumAccessibleSurfaceAreaDelta;

    /*
     description of the original residue + environment without any comparison on mutated site
     */
    private final double anyMutationFrequency;
    private final boolean anyMutationObserved;
    private final double deletionFrequency;
    private final boolean deletionObserved;
    private final double mutationToSimilarGutteridgeFrequency;
    private final boolean mutationToSimilarGutteridgeObserved;
    private final double mutationToTargetFrequency;
    private final boolean mutationToTargetObserved;
    private final double bindingSiteFrequency;
    private final boolean bindingSiteObserved;
    private final double siteOfInterestFrequency;
    private final boolean siteOfInterestObserved;
    private final double mutagenFrequency;
    private final boolean mutagenObserved;
    private final double variantFrequency;
    private final boolean variantObserved;
    private final double disulfidFrequency;
    private final boolean disulfidObserved;
    private final boolean disulfidCapabilitiesLost;

    /*
     environment delta
     */
    private final FeatureVector environmentDelta;
    /*
     amino acid delta
     */
    private final FeatureVector aminoAcidDelta;

    public MutationDescriptorImpl(String sequenceIdentifier,
                                  String mutationDescription,
                                  AminoAcid originalAminoAcid,
                                  AminoAcid mutatedAminoAcid,
                                  FeatureVector environmentDelta,
                                  FeatureVector aminoAcidDelta,
                                  double sequenceCount,
                                  double structureCount) {
        this.sequenceIdentifier = sequenceIdentifier;
        this.mutationDescription = mutationDescription;

        this.originalGutteridge = originalAminoAcid.getGroupPrototype().getGutteridgeGrouping();
        this.mutatedGutteridge = mutatedAminoAcid.getGroupPrototype().getGutteridgeGrouping();
        this.originalAminoAcidFunctional = !originalGutteridge.equals(GroupPrototype.GutteridgeGrouping.NONE);
        this.mutatedAminoAcidFunctional = !mutatedGutteridge.equals(GroupPrototype.GutteridgeGrouping.NONE);
        this.mutationChangedGutteridge = originalGutteridge != mutatedGutteridge;

        this.maximumAccessibleSurfaceAreaDelta = originalAminoAcid.getGroupPrototype().getMaximumAccessibleSurfaceArea() - mutatedAminoAcid.getGroupPrototype().getMaximumAccessibleSurfaceArea();

        AminoAcid.Family originalAminoAcidFamily = AminoAcid.Family.resolveGroupPrototype(originalAminoAcid.getGroupPrototype());
        AminoAcid.Family mutatedAminoAcidFamily = AminoAcid.Family.resolveGroupPrototype(mutatedAminoAcid.getGroupPrototype());
        FeatureContainer originalAminoAcidFeatureContainer = originalAminoAcid.getFeatureContainer();
        SequenceConservationProfile originalAminoAcidSequenceConservationProfile = originalAminoAcidFeatureContainer.getFeature(SequenceConservationProfile.class);
        this.anyMutationFrequency = originalAminoAcidSequenceConservationProfile.getFrequenceTable()
                .entrySet()
                .stream()
                .filter(entry -> entry.getKey() != originalAminoAcidFamily)
                .mapToDouble(Map.Entry::getValue)
                .sum();
        this.anyMutationObserved = anyMutationFrequency > 0;
        this.deletionFrequency = originalAminoAcidSequenceConservationProfile.getDeletionFrequency();
        this.deletionObserved = deletionFrequency > 0;
        this.mutationToSimilarGutteridgeFrequency = originalAminoAcidSequenceConservationProfile.getFrequenceTable()
                .entrySet()
                .stream()
                .filter(entry -> entry.getKey().getGutteridgeGrouping() == originalGutteridge)
                .mapToDouble(Map.Entry::getValue)
                .sum();
        this.mutationToSimilarGutteridgeObserved = mutationToSimilarGutteridgeFrequency > 0;
        this.mutationToTargetFrequency = originalAminoAcidSequenceConservationProfile.getFrequenceTable()
                .entrySet()
                .stream()
                .filter(entry -> entry.getKey() == originalAminoAcidFamily)
                .mapToDouble(Map.Entry::getValue)
                .sum();
        this.mutationToTargetObserved = mutationToTargetFrequency > 0;

        UniProtFeatureContainer featureContainer = originalAminoAcid.getFeatureContainer().getFeature(UniProtFeatureContainer.class);
        this.bindingSiteFrequency = featureContainer.getFeatures(FeatureType.CA_BIND,
                FeatureType.ZN_FING,
                FeatureType.DNA_BIND,
                FeatureType.NP_BIND,
                FeatureType.METAL,
                FeatureType.BINDING,
                FeatureType.LIPID)
                .size() / sequenceCount;
        this.bindingSiteObserved = bindingSiteFrequency > 0;
        this.siteOfInterestFrequency = featureContainer.getFeatures(FeatureType.REGION,
                FeatureType.ACT_SITE,
                FeatureType.SITE)
                .size() / sequenceCount;
        this.siteOfInterestObserved = siteOfInterestFrequency > 0;
        this.mutagenFrequency = featureContainer.getFeatures(FeatureType.MUTAGEN).size() / sequenceCount;
        this.mutagenObserved = mutagenFrequency > 0;
        this.variantFrequency = featureContainer.getFeatures(FeatureType.VAR_SEQ,
                FeatureType.VARIANT).size() / sequenceCount;
        this.variantObserved = variantFrequency > 0;
        this.disulfidFrequency = featureContainer.getFeatures(FeatureType.DISULFID).size() / sequenceCount;
        this.disulfidObserved = disulfidFrequency > 0;
        this.disulfidCapabilitiesLost = originalAminoAcidFamily == AminoAcid.Family.CYSTEINE && mutatedAminoAcidFamily != AminoAcid.Family.CYSTEINE;

        this.environmentDelta = environmentDelta;
        this.aminoAcidDelta = aminoAcidDelta;
    }

    public String getSequenceIdentifier() {
        return sequenceIdentifier;
    }

    public String getMutationDescription() {
        return mutationDescription;
    }

    public GroupPrototype.GutteridgeGrouping getOriginalGutteridge() {
        return originalGutteridge;
    }

    public GroupPrototype.GutteridgeGrouping getMutatedGutteridge() {
        return mutatedGutteridge;
    }

    public boolean isOriginalAminoAcidFunctional() {
        return originalAminoAcidFunctional;
    }

    public boolean isMutatedAminoAcidFunctional() {
        return mutatedAminoAcidFunctional;
    }

    public boolean isMutationChangedGutteridge() {
        return mutationChangedGutteridge;
    }

    public double getMaximumAccessibleSurfaceAreaDelta() {
        return maximumAccessibleSurfaceAreaDelta;
    }

    public boolean isAnyMutationObserved() {
        return anyMutationObserved;
    }

    public double getAnyMutationFrequency() {
        return anyMutationFrequency;
    }

    public boolean isDeletionObserved() {
        return deletionObserved;
    }

    public double getDeletionFrequency() {
        return deletionFrequency;
    }

    public boolean isMutationToSimilarGutteridgeObserved() {
        return mutationToSimilarGutteridgeObserved;
    }

    public double getMutationToSimilarGutteridgeFrequency() {
        return mutationToSimilarGutteridgeFrequency;
    }

    public boolean isMutationToTargetObserved() {
        return mutationToTargetObserved;
    }

    public double getMutationToTargetFrequency() {
        return mutationToTargetFrequency;
    }

    public double getBindingSiteFrequency() {
        return bindingSiteFrequency;
    }

    public boolean isBindingSiteObserved() {
        return bindingSiteObserved;
    }

    public double getSiteOfInterestFrequency() {
        return siteOfInterestFrequency;
    }

    public boolean isSiteOfInterestObserved() {
        return siteOfInterestObserved;
    }

    public double getMutagenFrequency() {
        return mutagenFrequency;
    }

    public boolean isMutagenObserved() {
        return mutagenObserved;
    }

    public double getVariantFrequency() {
        return variantFrequency;
    }

    public boolean isVariantObserved() {
        return variantObserved;
    }

    public double getDisulfidFrequency() {
        return disulfidFrequency;
    }

    public boolean isDisulfidObserved() {
        return disulfidObserved;
    }

    public boolean isDisulfidCapabilitiesLost() {
        return disulfidCapabilitiesLost;
    }

    public FeatureVector getEnvironmentDelta() {
        return environmentDelta;
    }

    public double getEnvironmentDeltaRasa() {
        return environmentDelta.getRasa();
    }

    public double getEnvironmentDeltaEnergy() {
        return environmentDelta.getEnergy();
    }

    public double getEnvironmentDeltaLigandContacts() {
        return environmentDelta.getLigandContacts();
    }

    public double getEnvironmentDeltaLoopFraction() {
        return environmentDelta.getLoopFraction();
    }

    public FeatureVector getAminoAcidDelta() {
        return aminoAcidDelta;
    }

    public double getAminoAcidDeltaRasa() {
        return aminoAcidDelta.getRasa();
    }

    public double getAminoAcidDeltaEnergy() {
        return aminoAcidDelta.getEnergy();
    }

    public double getAminoAcidDeltaLigandContacts() {
        return aminoAcidDelta.getLigandContacts();
    }

    public double getAminoAcidDeltaLoopFraction() {
        return aminoAcidDelta.getLoopFraction();
    }

    public double[] asDoubleVector() {
        return Stream.of(originalAminoAcidFunctional,
                mutatedAminoAcidFunctional,
                mutationChangedGutteridge,
                maximumAccessibleSurfaceAreaDelta,
                anyMutationFrequency,
                anyMutationObserved,
                deletionFrequency,
                deletionObserved,
                mutationToSimilarGutteridgeFrequency,
                mutationToSimilarGutteridgeObserved,
                mutationToTargetFrequency,
                mutationToTargetObserved,
                bindingSiteFrequency,
                bindingSiteObserved,
                siteOfInterestFrequency,
                siteOfInterestObserved,
                mutagenFrequency,
                mutagenObserved,
                variantFrequency,
                variantObserved,
                disulfidFrequency,
                disulfidObserved,
                disulfidCapabilitiesLost,
                environmentDelta.getRasa(),
                environmentDelta.getEnergy(),
                environmentDelta.getLigandContacts(),
                environmentDelta.getLoopFraction(),
                aminoAcidDelta.getRasa(),
                aminoAcidDelta.getEnergy(),
                aminoAcidDelta.getLigandContacts(),
                aminoAcidDelta.getLoopFraction())
                .mapToDouble(this::mapToDouble)
                .toArray();
    }

    private double mapToDouble(Object object) {
        if(object instanceof Boolean) {
            return (Boolean) object ? 1.0 : 0.0;
        }
        if(object instanceof Number) {
            return ((Number) object).doubleValue();
        }
        throw new IllegalArgumentException("object is of unexpected class '" + object.getClass().getSimpleName() + "'");
    }

    @Override
    public String asArffString() {
        return sequenceIdentifier + "," +
                mutationDescription + "," +
                DoubleStream.of(asDoubleVector())
                .mapToObj(StandardNumberFormat::format)
                .collect(Collectors.joining(","));
    }

    @Override
    public String toString() {
        return asArffString();
    }
}
