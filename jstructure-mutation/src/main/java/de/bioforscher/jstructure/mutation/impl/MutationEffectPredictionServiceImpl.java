package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.impl.ClustalOmegaRestQuery;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologyAnnotator;
import de.bioforscher.jstructure.mmm.MacromolecularMinerBridge;
import de.bioforscher.jstructure.mmm.impl.MacromolecularMinerBridgeImpl;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.mutation.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.kraken.interfaces.common.Value;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The implementation of the mutation effect prediction service.
 * Created by bittrich on 7/14/17.
 */
public class MutationEffectPredictionServiceImpl implements MutationEffectPredictionService {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionServiceImpl.class);
    private final UniProtHomologyAnnotator uniProtHomologyAnnotator;
    private final MultipleSequenceAligner multipleSequenceAligner;
    private final LigandContactScreener ligandContactScreener;
    private final MutatorService mutatorService;
    private final SequenceConservationCalculator sequenceConservationCalculator;
    private final MacromolecularMinerBridge macromolecularMinerBridge;
    private final AbstractFeatureProvider accessibleSurfaceCalculator;
    private final AbstractFeatureProvider energyProfileCalculator;
    private final AbstractFeatureProvider loopFractionCalculator;

    public MutationEffectPredictionServiceImpl() {
        this(new ClustalOmegaRestQuery(),
                new LigandContactScreenerImpl(),
                new ScwrlMutatorServiceImpl(),
                new SequenceConservationCalculatorImpl(),
                new MacromolecularMinerBridgeImpl());
    }

    public MutationEffectPredictionServiceImpl(MultipleSequenceAligner multipleSequenceAligner,
                                               LigandContactScreener ligandContactScreener,
                                               MutatorService mutatorService,
                                               SequenceConservationCalculator sequenceConservationCalculator,
                                               MacromolecularMinerBridge macromolecularMinerBridge) {
        this.multipleSequenceAligner = multipleSequenceAligner;
        this.ligandContactScreener = ligandContactScreener;
        this.mutatorService = mutatorService;
        this.sequenceConservationCalculator = sequenceConservationCalculator;
        this.macromolecularMinerBridge = macromolecularMinerBridge;
        this.uniProtHomologyAnnotator = new UniProtHomologyAnnotator();
        this.accessibleSurfaceCalculator = FeatureProviderRegistry.resolve(AccessibleSurfaceArea.class);
        this.energyProfileCalculator = FeatureProviderRegistry.resolve(EnergyProfile.class);
        this.loopFractionCalculator = FeatureProviderRegistry.resolve(LoopFraction.class);
    }

    @Override
    public MutationJob createMutationJob(String jobName, Chain referenceChain) {
        MutationJobImpl mutationJob = new MutationJobImpl(jobName, referenceChain);
        String uuid = mutationJob.getUuid().toString();
        logger.info("[{}] started job '{}' with reference chain '{}' at {}",
                uuid,
                mutationJob.getJobName(),
                mutationJob.getReferenceChain().getChainIdentifier().getFullName(),
                mutationJob.getCreationTime());

        // compute features of reference chain, track exceptions to enable fail-fast behaviour
        logger.info("[{}] computing features for reference chain",
                mutationJob.getUuid());
        List<Exception> exceptionsDuringReferenceChainAnnotation = annotateChainAndReportExceptions(referenceChain);
        if(exceptionsDuringReferenceChainAnnotation.size() > 0) {
            exceptionsDuringReferenceChainAnnotation
                    .forEach(ex -> logger.warn("[{}] encountered exception during annotation of refernce chain",
                            mutationJob.getUuid(),
                            ex));
            //TODO decide to what degree exceptions here will compromise the job
        }

        // execute BLAST query against SWISS-PROT
        String sequence = mutationJob.getReferenceChain().getAminoAcidSequence();
        logger.info("[{}] executing BLAST query with sequence {}",
                uuid,
                sequence);
        uniProtHomologyAnnotator.process(mutationJob.getReferenceChain().getParentProtein());

        // fetch homologous UniProt entries
        UniProtHomologousEntryContainer uniProtHomologousEntryContainer = mutationJob.getReferenceChain().getFeatureContainer().getFeature(UniProtHomologousEntryContainer.class);
        mutationJob.setHomologousSequences(uniProtHomologousEntryContainer.getUniProtEntries());
        logger.info("[{}] homologous sequences {}",
                uuid,
                uniProtHomologousEntryContainer.getUniProtEntries()
                .stream()
                .map(UniProtEntry::getPrimaryUniProtAccession)
                .map(Value::getValue)
                .collect(Collectors.toList()));

        // fetch homologous PDB entries
        List<Chain> homologousPdbChains = uniProtHomologousEntryContainer.getHomologousChains()
                .stream()
                // ignore the reference chain itself
                .filter(chainIdentifier -> !chainIdentifier.equals(mutationJob.getReferenceChain().getChainIdentifier()))
                .map(chainIdentifier -> fetchProteinChain(chainIdentifier, uuid))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());

        // compute feature for homologous chains - remove compromised entries
        List<Chain> homologousChainsWhoseAnnotationFailed = homologousPdbChains.stream()
                .filter(chain -> annotateChainAndReportExceptions(chain).size() > 0)
                .collect(Collectors.toList());
        homologousChainsWhoseAnnotationFailed.forEach(chain -> logger.warn("[{}] feature annotation failed for homologous chain {}",
                mutationJob.getUuid(),
                chain.getChainIdentifier().getFullName()));
        homologousPdbChains.removeAll(homologousChainsWhoseAnnotationFailed);
        logger.info("[{}] homologous PDB chains {}",
                uuid,
                homologousPdbChains
                        .stream()
                        .map(Chain::getChainIdentifier)
                        .map(ChainIdentifier::getFullName)
                        .collect(Collectors.toList()));
        mutationJob.setHomologousPdbChains(homologousPdbChains);

        // create multiple-sequence alignment of everything
        renumberSequences(mutationJob);

        return mutationJob;
    }

    private void renumberSequences(MutationJob mutationJob) {
        // execute ClustalOmega on all sequences
        List<String> sequences = new ArrayList<>();
        Chain referenceChain = mutationJob.getReferenceChain();

        sequences.add(">" + referenceChain.getChainIdentifier().getFullName() + System.lineSeparator() + referenceChain.getAminoAcidSequence());
        mutationJob.getHomologousSequences()
                .stream()
                .map(entry -> ">" + entry.getPrimaryUniProtAccession().getValue() + System.lineSeparator() + entry.getSequence().getValue())
                .forEach(sequences::add);
        mutationJob.getHomologousPdbChains()
                .stream()
                .map(chain -> ">" + chain.getChainIdentifier().getFullName() + System.lineSeparator() + chain.getAminoAcidSequence())
                .forEach(sequences::add);

        logger.info("[{}] executing multiple-sequence alignment on {} sequence by {}",
                mutationJob.getUuid(),
                sequences.size(),
                multipleSequenceAligner.getClass().getSimpleName());
        Map<String, String> alignmentMap = multipleSequenceAligner.align(sequences).getAlignedSequences();
        logger.info("[{}] renumbering sequences",
                mutationJob.getUuid());
        getStreamOfAllChains(mutationJob)
                .forEach(chain -> {
                    String alignmentString = alignmentMap.get(chain.getChainIdentifier().getFullName());
                    int alignmentLength = alignmentString.length();
                    int consumedAminoAcidsInChain = 0;
                    List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
                    for(int sequencePosition = 0; sequencePosition < alignmentLength; sequencePosition++) {
                        try {
                            char characterInAlignment = alignmentString.charAt(sequencePosition);
                            if (characterInAlignment == '-') {
                                continue;
                            }
                            AminoAcid aminoAcid = aminoAcids.get(consumedAminoAcidsInChain);
                            consumedAminoAcidsInChain++;
                            // renumber to position in alignment string (+1 for classic Java offset)
                            aminoAcid.setResidueIdentifier(IdentifierFactory.createResidueIdentifier(sequencePosition + 1));
                        } catch (Exception e) {
                            // missing amino acids or premature end of List
                        }
                    }
                });
    }

    /**
     * Tries to annotate a chain given a defined set of features. The behaviour is that all providers will get their try
     * in processing the protein. Some may fail: if so, results of others should be still sane, but the function will
     * not an empty list but all exceptions caught.
     * @param chain the (isolated) chain to process
     * @return a collection of exceptions which arose during computation
     */
    private List<Exception> annotateChainAndReportExceptions(Chain chain) {
        return Stream.of(accessibleSurfaceCalculator,
                energyProfileCalculator,
                loopFractionCalculator)
                .map(featureProvider -> annotateChainAndReportPotentialException(chain, featureProvider))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }

    /**
     * A dedicated way to run a feature provider on a single chain.
     * @param chain the chain to process
     * @param featureProvider the feature provider to employ
     * @return an optional containing an exception if one arose during computation
     */
    private Optional<Exception> annotateChainAndReportPotentialException(Chain chain, AbstractFeatureProvider featureProvider) {
        try {
            featureProvider.process(chain.getParentProtein());
            return Optional.empty();
        } catch (Exception e) {
            return Optional.of(e);
        }
    }

    /**
     * Convenience function to get access to all chains of a container.
     * @param mutationJob the container
     * @return all {@link Chain} instances associated with it - essentially just the reference and all homologous concat
     */
    private Stream<Chain> getStreamOfAllChains(MutationJob mutationJob) {
        return Stream.concat(Stream.of(mutationJob.getReferenceChain()), mutationJob.getHomologousPdbChains().stream());
    }

    /**
     * Internal class to handle the renumbering without compromising the original numbering.
     */
    public static class RenumberedResidueNumber extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Integer> {
        private final int renumberedResidueNumber;

        public RenumberedResidueNumber(int renumberedResidueNumber) {
            super(null);
            this.renumberedResidueNumber = renumberedResidueNumber;
        }

        public int getRenumberedResidueNumber() {
            return renumberedResidueNumber;
        }

        @Override
        public Integer getValue() {
            return renumberedResidueNumber;
        }

        /**
         * Convenience function to retrieve a particular amino acid by its renumbered identifier.
         * @param chain the chain to search in
         * @param renumberedResidueNumber the renumbered value the amino acid is supposed to feature
         * @return a optional containing a value when a suitable amino acid was found
         */
        public static Optional<AminoAcid> getAminoAcidByRenumberedResidueNumber(Chain chain, int renumberedResidueNumber) {
            return chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeatureContainer()
                            .getFeatureOptional(RenumberedResidueNumber.class)
                            .map(RenumberedResidueNumber::getRenumberedResidueNumber)
                            .map(rrr -> rrr == renumberedResidueNumber)
                            // return false when operation failed
                            .orElse(false))
                    .findFirst();
        }
    }

    private Optional<Chain> fetchProteinChain(ChainIdentifier chainIdentifier, String uuid) {
        try {
            Protein originalProtein = ProteinParser.source(chainIdentifier.getProteinIdentifier().getPdbId())
                    .minimalParsing(true)
                    .parse();
            Chain isolatedChain = originalProtein.select()
                    .chainName(chainIdentifier.getChainId())
                    .asChain();
            Protein isolatedProtein = new Protein(isolatedChain.getChainIdentifier().getProteinIdentifier());
            // will set the parent reference to the isolated chain
            isolatedProtein.addChain(isolatedChain);
            return Optional.of(isolatedChain);
        } catch (Exception e) {
            logger.info("[{}] failed to fetch protein '{}'",
                    uuid,
                    chainIdentifier.getFullName(),
                    e);
            return Optional.empty();
        }
    }

    @Override
    public MutationFeatureVector createMutationFeatureVector(MutationJob mutationJob, ResidueIdentifier residueIdentifierToMutate, AminoAcid.Family mutationTarget) {
        return null;
    }
}
