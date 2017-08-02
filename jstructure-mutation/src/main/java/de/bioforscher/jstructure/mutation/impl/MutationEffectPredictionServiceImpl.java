package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.align.MultipleSequenceAligner;
import de.bioforscher.jstructure.align.impl.ClustalOmegaWrapper;
import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformationCalculator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.MutationEffectPredictionService;
import de.bioforscher.jstructure.mutation.MutationFeatureVector;
import de.bioforscher.jstructure.mutation.MutationJob;
import de.bioforscher.jstructure.mutation.MutatorService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The implementation of the mutation effect prediction service.
 * Created by bittrich on 7/14/17.
 */
public class MutationEffectPredictionServiceImpl implements MutationEffectPredictionService {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionServiceImpl.class);
    private static final Pattern CHAIN_PATTERN = Pattern.compile("/");
    /**
     * The size of the sphere extracted around an amino acid to describe its surroundings.
     */
    private static final double AMINO_ACID_ENVIRONMENT_SIZE = 8.0;

    private final LocalBlastWrapper localBlastWrapper;
    private final EvolutionaryInformationCalculator evolutionaryInformationCalculator;
    private final MultipleSequenceAligner multipleSequenceAligner;
    private final MutatorService mutatorService;
    private final FeatureProvider accessibleSurfaceCalculator;
    private final FeatureProvider energyProfileCalculator;
    private final FeatureProvider loopFractionCalculator;

    public MutationEffectPredictionServiceImpl() {
        this(new LocalBlastWrapper(),
                new EvolutionaryInformationCalculator(),
                new ClustalOmegaWrapper(),
                new ScwrlMutatorServiceImpl());
    }

    public MutationEffectPredictionServiceImpl(LocalBlastWrapper localBlastWrapper,
                                               EvolutionaryInformationCalculator evolutionaryInformationCalculator,
                                               MultipleSequenceAligner multipleSequenceAligner,
                                               MutatorService mutatorService) {
        this.localBlastWrapper = localBlastWrapper;
        this.evolutionaryInformationCalculator = evolutionaryInformationCalculator;
        this.multipleSequenceAligner = multipleSequenceAligner;
        this.mutatorService = mutatorService;
        this.accessibleSurfaceCalculator = new AccessibleSurfaceAreaCalculator();
        this.energyProfileCalculator = new EnergyProfileCalculator();
        this.loopFractionCalculator = new LoopFractionCalculator();
    }

    @Override
    public MutationJob createMutationJob(String jobName, Chain referenceChain, LocalBlastWrapper.PsiBlastResult uniref50result) {
        MutationJobImpl mutationJob = new MutationJobImpl(jobName, referenceChain);
        String uuid = mutationJob.getUuid().toString();
        logger.info("[{}] started job '{}' with reference chain '{}' at {}",
                mutationJob.getJobName(),
                uuid,
                mutationJob.getReferenceChain().getChainIdentifier().getFullName(),
                mutationJob.getCreationTime());

        // compute features of reference chain, track exceptions to enable fail-fast behaviour
        logger.info("[{}] computing features for reference chain",
                mutationJob.getJobName());
        List<Exception> exceptionsDuringReferenceChainAnnotation = annotateChainAndReportExceptions(mutationJob.getReferenceChain());
        if(exceptionsDuringReferenceChainAnnotation.size() > 0) {
            exceptionsDuringReferenceChainAnnotation
                    .forEach(ex -> logger.warn("[{}] encountered exception during annotation of reference chain",
                            mutationJob.getUuid(),
                            ex));
            //TODO decide to what degree exceptions here will compromise the job
        }

        // assign PSSM conservation score to each amino acid
        evolutionaryInformationCalculator.assignPsiBlastResultToChain(referenceChain, uniref50result);

//        // execute PSI-BLAST query against SWISS-PROT
//        String sequence = mutationJob.getReferenceChain().getAminoAcidSequence();
//        logger.info("[{}] executing PSI-BLAST query against SWISS-PROT with sequence {}",
//                mutationJob.getJobName(),
//                sequence);
//        LocalBlastWrapper.PsiBlastResult swissprotResult = localBlastWrapper.executePsiBlastSwissProt(">query" + System.lineSeparator() + sequence);
//
//        List<UniProtEntry> uniProtEntries = swissprotResult.getAccessions()
//                .stream()
//                .map(this::getUniProtEntry)
//                .filter(Optional::isPresent)
//                .map(Optional::get)
//                .collect(Collectors.toList());
//        logger.info("[{}] homologous sequences {}",
//                mutationJob.getJobName(),
//                uniProtEntries.stream()
//                .map(UniProtEntry::getPrimaryUniProtAccession)
//                .map(Value::getValue)
//                .collect(Collectors.toList()));
//        mutationJob.setHomologousSequences(uniProtEntries);
//
//        // validate number of homologous sequences
//        if(uniProtEntries.isEmpty()) {
//            throw new IllegalArgumentException("could not find any sequence homologs");
//        }
//
//        List<ChainIdentifier> homologousChains = uniProtEntries.stream()
//                // take all referenced chains
//                .map(uniProtEntry -> uniProtEntry.getDatabaseCrossReferences(DatabaseType.PDB))
//                .flatMap(Collection::stream)
//                .flatMap(databaseCrossReference -> {
//                    //TODO could use range which also available in the format fourth -> "A=321-419"
//                    String[] split = databaseCrossReference.getFourth().getValue().split("=");
//                    // some entries reference multiple chains separated by '/'
//                    return CHAIN_PATTERN.splitAsStream(split[0])
//                            .map(chainId -> IdentifierFactory.createChainIdentifier(databaseCrossReference.getPrimaryId().getValue(), chainId));
//                })
//                .distinct()
//                .collect(Collectors.toList());
//        logger.info("[{}] homologous structure candidates: {}",
//                uuid,
//                homologousChains);
//        homologousChains = removeRedundancy(homologousChains);
//
//        // fetch homologous PDB entries
//        List<Chain> homologousPdbChains = homologousChains.stream()
//                // ignore the reference chain itself
//                .filter(chainIdentifier -> !chainIdentifier.equals(mutationJob.getReferenceChain().getChainIdentifier()))
//                .map(chainIdentifier -> fetchProteinChain(chainIdentifier, mutationJob.getJobName()))
//                .filter(Optional::isPresent)
//                .map(Optional::get)
//                .collect(Collectors.toList());
//
//        // compute feature for homologous chains - remove compromised entries
//        logger.info("[{}] computing features for homologous chains",
//                mutationJob.getJobName());
//        List<Chain> homologousChainsWhoseAnnotationFailed = homologousPdbChains.stream()
//                .filter(chain -> annotateChainAndReportExceptions(chain).size() > 0)
//                .collect(Collectors.toList());
//        homologousChainsWhoseAnnotationFailed.forEach(chain -> logger.warn("[{}] feature annotation failed for homologous chain {}",
//                mutationJob.getJobName(),
//                chain.getChainIdentifier().getFullName()));
//        homologousPdbChains.removeAll(homologousChainsWhoseAnnotationFailed);
//        logger.info("[{}] homologous PDB chains {}",
//                mutationJob.getJobName(),
//                homologousPdbChains
//                        .stream()
//                        .map(Chain::getChainIdentifier)
//                        .map(ChainIdentifier::getFullName)
//                        .collect(Collectors.toList()));
//        mutationJob.setHomologousPdbChains(homologousPdbChains);
//
//        // create multiple-sequence alignment of everything
//        renumberSequences(mutationJob);

        return mutationJob;
    }

//    private List<ChainIdentifier> removeRedundancy(List<ChainIdentifier> homologousChains) {
//        Map<ChainIdentifier, ChainIdentifier> processedChainIdentifiers = new HashMap<>();
//        List<ChainIdentifier> selectedChainIdentifiers = new ArrayList<>();
//
//        for(ChainIdentifier chainIdentifier : homologousChains) {
//            // already processed
//            if(processedChainIdentifiers.containsKey(chainIdentifier)) {
//                ChainIdentifier selectedChainIdentifier = processedChainIdentifiers.get(chainIdentifier);
//                if(selectedChainIdentifier.equals(chainIdentifier)) {
//                    selectedChainIdentifiers.add(chainIdentifier);
//                }
//                continue;
//            }
//
//            String url = "https://www.rcsb.org/pdb/rest/sequenceCluster?cluster=95&structureId=" +
//                            chainIdentifier.getProteinIdentifier().getPdbId() + "." + chainIdentifier.getChainId();
//            try {
//                Elements cluster = Jsoup.connect(url)
//                        .get()
//                        .getElementsByTag("pdbChain");
//                String[] representativeSplit = cluster
//                        .first()
//                        .attr("name")
//                        .split("\\.");
//                ChainIdentifier representative = IdentifierFactory.createChainIdentifier(representativeSplit[0],
//                        representativeSplit[1]);
//                for (Element element : cluster) {
//                    String[] elementSplit = element.attr("name")
//                            .split("\\.");
//                    processedChainIdentifiers.put(IdentifierFactory.createChainIdentifier(elementSplit[0],
//                            elementSplit[1]), representative);
//                }
//
//                if(representative.equals(chainIdentifier)) {
//                    selectedChainIdentifiers.add(chainIdentifier);
//                }
//            } catch (Exception e) {
////                logger.warn("[] failed to fetch sequence cluster for {} from pdb",
////                        chainIdentifier);
//            }
//        }
//        if(selectedChainIdentifiers.isEmpty()) {
//            return homologousChains;
//        }
//        return selectedChainIdentifiers;
//    }
//
//    private Optional<UniProtEntry> getUniProtEntry(String accession) {
//        try {
//            Query query = UniProtQueryBuilder.accession(accession);
//            QueryResult<UniProtEntry> result = uniProtService.getEntries(query);
//            return Optional.of(result.getFirstResult());
//        } catch (Exception e) {
//            logger.warn("could not retrieve UniProt entry for {}",
//                    accession,
//                    e);
//            return Optional.empty();
//        }
//    }
//
//    private void renumberSequences(MutationJobImpl mutationJob) {
//        // execute ClustalOmega on all sequences
//        List<String> sequences = new ArrayList<>();
//        Chain referenceChain = mutationJob.getReferenceChain();
//
//        sequences.add(">" + referenceChain.getChainIdentifier().getFullName() + System.lineSeparator() + referenceChain.getAminoAcidSequence());
//        mutationJob.getHomologousSequences()
//                .stream()
//                .map(entry -> ">" + entry.getPrimaryUniProtAccession().getValue() + System.lineSeparator() + entry.getSequence().getValue())
//                .forEach(sequences::add);
//        mutationJob.getHomologousPdbChains()
//                .stream()
//                .map(chain -> ">" + chain.getChainIdentifier().getFullName() + System.lineSeparator() + chain.getAminoAcidSequence())
//                .forEach(sequences::add);
//
//        logger.info("[{}] executing multiple-sequence alignment on {} sequence by {}",
//                mutationJob.getJobName(),
//                sequences.size(),
//                multipleSequenceAligner.getClass().getSimpleName());
//        Map<String, String> alignmentMap = multipleSequenceAligner.align(sequences).getAlignedSequences();
//        logger.info("[{}] renumbering sequences",
//                mutationJob.getJobName());
//        getStreamOfAllChains(mutationJob)
//                .forEach(chain -> {
//                    String alignmentString = alignmentMap.get(chain.getChainIdentifier().getFullName());
//                    int alignmentLength = alignmentString.length();
//                    int consumedAminoAcidsInChain = 0;
//                    List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
//                    for(int sequencePosition = 0; sequencePosition < alignmentLength; sequencePosition++) {
//                        try {
//                            char characterInAlignment = alignmentString.charAt(sequencePosition);
//                            if (characterInAlignment == '-') {
//                                continue;
//                            }
//                            AminoAcid aminoAcid = aminoAcids.get(consumedAminoAcidsInChain);
//                            consumedAminoAcidsInChain++;
//                            // renumber to position in alignment string (+1 for classic Java offset)
//                            aminoAcid.getFeatureContainer().addFeature(new RenumberedResidueNumber(sequencePosition + 1));
//                        } catch (Exception e) {
//                            // missing amino acids or premature end of List
//                        }
//                    }
//                });
//
//        mutationJob.setAlignmentMap(alignmentMap);
//
//        // map UniProt features onto reference sequence
//        annotateReferenceChainWithUniProtFeatures(referenceChain, mutationJob);
//    }
//
//    private void annotateReferenceChainWithUniProtFeatures(Chain referenceChain, MutationJob mutationJob) {
//        Map<String, String> alignmentMap = mutationJob.getAlignmentMap();
//        String referenceChainAlignmentString = alignmentMap.get(referenceChain.getChainIdentifier().getFullName());
//        List<AminoAcid> aminoAcids = referenceChain.aminoAcids().collect(Collectors.toList());
//        // initialize all amino acids as empty annotations
//        aminoAcids.forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new UniProtFeatureContainer(null)));
//
//        mutationJob.getHomologousSequences()
//                .forEach(uniProtEntry -> {
//                    String accession = uniProtEntry.getPrimaryUniProtAccession().getValue();
//                    String uniProtEntryAlignmentString = alignmentMap.get(accession);
//
//                    for(int alignmentStringIndex = 0; alignmentStringIndex < referenceChainAlignmentString.length(); alignmentStringIndex++) {
//                        char referenceChar = referenceChainAlignmentString.charAt(alignmentStringIndex);
//                        // continue if gap in reference sequence
//                        if(referenceChar == '-') {
//                            continue;
//                        }
//
//                        String matchSequence = uniProtEntryAlignmentString.substring(0, alignmentStringIndex);
//                        int matchGaps = countGaps(matchSequence);
//                        int positionInMatchSequence = alignmentStringIndex + 1 - matchGaps;
//
//                        List<Feature> features = getFeaturesForCurrentPosition(uniProtEntry, positionInMatchSequence);
//                        // continue if no feature at current position of match
//                        if(features.isEmpty()) {
//                            continue;
//                        }
//
//                        // operate on partial string, ignore last position for now
//                        String referenceSequence = referenceChainAlignmentString.substring(0, alignmentStringIndex);
//                        int referenceGaps = countGaps(referenceSequence);
//                        int positionInReferenceSequence = alignmentStringIndex + 1 - referenceGaps;
//
//                        Optional<AminoAcid> aminoAcidOptional = getSafely(aminoAcids, positionInReferenceSequence - 1);
//                        if(aminoAcidOptional.isPresent()) {
//                            UniProtFeatureContainer featureContainer = aminoAcidOptional.get()
//                                    .getFeatureContainer()
//                                    .getFeature(UniProtFeatureContainer.class);
//                            features.forEach(feature -> featureContainer.addFeature(accession, feature));
//                        }
//                    }
//                });
//    }
//
//    private Optional<AminoAcid> getSafely(List<AminoAcid> aminoAcids, int index) {
//        try {
//            return Optional.of(aminoAcids.get(index));
//        } catch (ArrayIndexOutOfBoundsException e) {
//            return Optional.empty();
//        }
//    }
//
//    private List<Feature> getFeaturesForCurrentPosition(UniProtEntry uniProtEntry, int positionInMatchSequence) {
//        return uniProtEntry.getFeatures()
//                .stream()
//                .filter(feature -> featureEmbedsCurrentPosition(feature, positionInMatchSequence))
//                .collect(Collectors.toList());
//    }
//
//    private boolean featureEmbedsCurrentPosition(Feature feature, int position) {
//        FeatureLocation featureLocation = feature.getFeatureLocation();
//        int start = featureLocation.getStart();
//        int end = (featureLocation.isEndAvailable() ? featureLocation.getEnd() : start);
//        return start <= position && position <= end;
//    }
//
//    private int countGaps(String sequence) {
//        int count = 0;
//        for(int i = 0; i < sequence.length(); i++) {
//            if(sequence.charAt(i) == '-') {
//                count++;
//            }
//        }
//        return count;
//    }

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
    private Optional<Exception> annotateChainAndReportPotentialException(Chain chain, FeatureProvider featureProvider) {
        try {
            featureProvider.process(chain.getParentStructure());
            return Optional.empty();
        } catch (Exception e) {
            return Optional.of(e);
        }
    }

//    /**
//     * Convenience function to get access to all chains of a container.
//     * @param mutationJob the container
//     * @return all {@link Chain} instances associated with it - essentially just the reference and all homologous concat
//     */
//    private Stream<Chain> getStreamOfAllChains(MutationJob mutationJob) {
//        return Stream.concat(Stream.of(mutationJob.getReferenceChain()), mutationJob.getHomologousPdbChains().stream());
//    }
//
//    /**
//     * Internal class to handle the renumbering without compromising the original numbering.
//     */
//    public static class RenumberedResidueNumber extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Integer> {
//        private final int renumberedResidueNumber;
//
//        public RenumberedResidueNumber(int renumberedResidueNumber) {
//            super(null);
//            this.renumberedResidueNumber = renumberedResidueNumber;
//        }
//
//        public int getRenumberedResidueNumber() {
//            return renumberedResidueNumber;
//        }
//
//        @Override
//        public Integer getValue() {
//            return renumberedResidueNumber;
//        }
//
//        /**
//         * Convenience function to retrieve a particular amino acid by its renumbered identifier.
//         * @param chain the chain to search in
//         * @param renumberedResidueNumber the renumbered value the amino acid is supposed to feature
//         * @return a optional containing a value when a suitable amino acid was found
//         */
//        public static Optional<AminoAcid> getAminoAcidByRenumberedResidueNumber(Chain chain, int renumberedResidueNumber) {
//            return chain.aminoAcids()
//                    .filter(aminoAcid -> aminoAcid.getFeatureContainer()
//                            .getFeatureOptional(RenumberedResidueNumber.class)
//                            .map(RenumberedResidueNumber::getRenumberedResidueNumber)
//                            .map(rrr -> rrr == renumberedResidueNumber)
//                            // return false when operation failed
//                            .orElse(false))
//                    .findFirst();
//        }
//    }
//
//    private Optional<Chain> fetchProteinChain(ChainIdentifier chainIdentifier, String jobName) {
//        try {
//            Structure originalProtein = StructureParser.source(chainIdentifier.getProteinIdentifier().getPdbId())
//                    .minimalParsing(true)
//                    .parse();
//            Chain isolatedChain = originalProtein.select()
//                    .chainName(chainIdentifier.getChainId())
//                    .asChain();
//            Structure isolatedProtein = new Structure(isolatedChain.getChainIdentifier().getProteinIdentifier());
//            // will set the parent reference to the isolated chain
//            isolatedProtein.addChain(isolatedChain);
//            return Optional.of(isolatedChain);
//        } catch (Exception e) {
//            logger.info("[{}] failed to fetch protein '{}'",
//                    jobName,
//                    chainIdentifier.getFullName(),
//                    e);
//            return Optional.empty();
//        }
//    }

    @Override
    public MutationFeatureVector createMutationFeatureVector(MutationJob mutationJob,
                                                             ChainIdentifier chainIdentifier,
                                                             ResidueIdentifier residueIdentifierToMutate,
                                                             AminoAcid.Family mutationTarget) {
        // original data
        Chain originalChain = mutationJob.getReferenceChain();
        AminoAcid originalAminoAcid = originalChain.select()
                .residueIdentifier(residueIdentifierToMutate)
                .asAminoAcid();

        // mutated data
        Structure mutatedProtein = mutatorService.mutateAminoAcid(mutationJob.getReferenceChain().getParentStructure(),
                chainIdentifier,
                residueIdentifierToMutate,
                mutationTarget);
        Chain mutatedChain = mutatedProtein.select()
                .chainIdentifier(chainIdentifier)
                .asChain();
        AminoAcid mutatedAminoAcid = mutatedChain.select()
                .residueIdentifier(residueIdentifierToMutate)
                .asAminoAcid();

        // compute features on mutated protein
        List<Exception> exceptions = annotateChainAndReportExceptions(mutatedChain);
        if(!exceptions.isEmpty()) {
            exceptions.forEach(ex -> logger.warn("[{}] exception during feature annotation of of mutated protein {}{}{}",
                    mutationJob.getJobName(),
                    originalAminoAcid.getOneLetterCode(),
                    residueIdentifierToMutate.getResidueNumber(),
                    mutatedAminoAcid.getOneLetterCode(),
                    ex));
        }

        logger.info("[{}] creating feature vector for {}: {}{}{}",
                mutationJob.getJobName(),
                chainIdentifier.getFullName(),
                originalAminoAcid.getOneLetterCode(),
                residueIdentifierToMutate.getResidueNumber(),
                mutatedAminoAcid.getOneLetterCode());

        // extract environments
        List<AminoAcid> originalEnvironment = extractEnvironment(originalAminoAcid, mutationJob.getJobName());
        List<AminoAcid> mutatedEnvironment = extractEnvironment(mutatedAminoAcid, mutationJob.getJobName());

        // compute feature vectors
        PhysicochemicalFeatureVector fvOriginalAminoAcid = new PhysicochemicalFeatureVector(originalAminoAcid);
        PhysicochemicalFeatureVector fvMutatedAminoAcid = new PhysicochemicalFeatureVector(mutatedAminoAcid);
        PhysicochemicalFeatureVector fvOriginalEnvironment = new PhysicochemicalFeatureVector(originalEnvironment);
        PhysicochemicalFeatureVector fvMutatedEnvironment = new PhysicochemicalFeatureVector(mutatedEnvironment);

        // compute delta
        DeltaPhysicochemicalFeatureVector dfvAminoAcid = new DeltaPhysicochemicalFeatureVector(fvOriginalAminoAcid, fvMutatedAminoAcid);
        DeltaPhysicochemicalFeatureVector dfvEnvironment = new DeltaPhysicochemicalFeatureVector(fvOriginalEnvironment, fvMutatedEnvironment);

        // general information
        return new MutationFeatureVectorImpl(mutationJob,
                chainIdentifier,
                residueIdentifierToMutate,

                // amino acid instances
                originalAminoAcid,
                mutatedAminoAcid,

                originalEnvironment,
                mutatedEnvironment,

                fvOriginalAminoAcid,

                // feature vectors
//                fvOriginalAminoAcid,
//                fvMutatedAminoAcid,
//                fvOriginalEnvironment,
//                fvMutatedEnvironment,

                // delta vectors
                dfvAminoAcid,
                dfvEnvironment
        );
    }

    /**
     * Extracts the environment (i.e. all amino acids of this chain at most 8 A away) and reports them as list.
     * @param aminoAcid the amino acid to investigate
     * @param jobName the job's name
     * @return a collection of all surrounding amino acids
     */
    private List<AminoAcid> extractEnvironment(AminoAcid aminoAcid, String jobName) {
        try {
            Atom referenceAtom = aminoAcid.getCa();
            return aminoAcid.getParentChain().select()
                    .aminoAcids()
                    .groupDistance(referenceAtom, AMINO_ACID_ENVIRONMENT_SIZE)
                    .asFilteredGroups()
                    .map(AminoAcid.class::cast)
                    .collect(Collectors.toList());
        } catch (Exception e) {
            logger.warn("could not extract environment around {}",
                    jobName,
                    aminoAcid,
                    e);
            return new ArrayList<>();
        }
    }

    /**
     * Describes the features of any number of amino acids.
     */
    static class PhysicochemicalFeatureVector {
        protected double rasa;
        protected double loopFraction;
        protected double energy;

        protected PhysicochemicalFeatureVector() {

        }

        PhysicochemicalFeatureVector(List<AminoAcid> aminoAcids) {
            this(aminoAcids.stream(),
                    aminoAcids.size(),
                    aminoAcids.get(0).getParentChain().getParentStructure());
        }

        PhysicochemicalFeatureVector(AminoAcid... aminoAcids) {
            this(Stream.of(aminoAcids),
                    aminoAcids.length,
                    aminoAcids[0].getParentChain().getParentStructure());
        }

        private PhysicochemicalFeatureVector(Stream<AminoAcid> aminoAcidStream,
                                             double numberOfAminoAcids,
                                             Structure structure) {
            aminoAcidStream
                    .forEach(aminoAcid -> {
                        rasa += aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea() / numberOfAminoAcids;
                        loopFraction += aminoAcid.getFeature(LoopFraction.class).getLoopFraction() / numberOfAminoAcids;
                        energy += aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy() / numberOfAminoAcids;
                    });
        }

        public double getRasa() {
            return rasa;
        }

        public double getLoopFraction() {
            return loopFraction;
        }

        public double getEnergy() {
            return energy;
        }
    }

    /**
     * Describes the change in features when a mutation occurs.
     */
    static class DeltaPhysicochemicalFeatureVector extends PhysicochemicalFeatureVector {
        DeltaPhysicochemicalFeatureVector(PhysicochemicalFeatureVector original, PhysicochemicalFeatureVector mutated) {
            this.rasa = original.rasa - mutated.rasa;
            this.loopFraction = original.loopFraction - mutated.loopFraction;
            this.energy = original.energy - mutated.energy;
        }
    }
}
