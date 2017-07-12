package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.mutation.MutationJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

/**
 * The implementation of the data structure.
 * Created by bittrich on 7/11/17.
 */
class MutationJobImpl implements MutationJob {
    private static final Logger logger = LoggerFactory.getLogger(MutationJobImpl.class);
    private static final double INTERACTION_CUTOFF = 8.0;
    private static final Path MOLECULAR_MINER_WORKING_PATH = Paths.get("/tmp/mmm/");

    /**
     * The name of this job.
     */
    private final String identifier;
    /**
     * The query sequence which was used to create this job.
     */
    private final String querySequence;
    /**
     * The query sequence wrapped as Protein and Chain object. Does not provide features or coordinates.
     */
    private final Protein queryProtein;
    private final Chain queryChain;
    /**
     * The PDB-chain chosen as reference. Does provide features and coordinates.
     */
    private Protein referenceProtein;
    private Chain referenceChain;
    /**
     * Access to the BLAST alignment and all hits and the data of the respective UniProt entries.
     */
    private UniProtHomologousEntryContainer homologousEntryContainer;
    /**
     * All similar structures in the PDB (likely to contain the reference structure again). Provide features.
     */
    private List<Chain> homologousPdbChains;
    /**
     * The alignment map (id, alignmentString) provides the MSA for the querySequence and all PDB chains.
     */
    private Map<String, String> alignmentMap;

    MutationJobImpl(String identifier, String querySequence) {
        this.identifier = identifier;
        this.querySequence = querySequence;
        this.queryProtein = createProtein(identifier, querySequence);
        this.queryChain = queryProtein.getChains().get(0);
    }

    /**
     * Creates a pseudo protein wrapping only the sequence.
     * @param identifier the name
     * @param sequence the sequence of groups to create
     * @return the hollow instance of a protein
     */
    static Protein createProtein(String identifier, String sequence) {
        ProteinIdentifier proteinIdentifier = IdentifierFactory.createProteinIdentifier("", identifier);
        Protein protein = new Protein(proteinIdentifier);
        Chain chain = new Chain(IdentifierFactory.createChainIdentifier(proteinIdentifier, "X"));
        for(int residueNumber = 1; residueNumber <= sequence.length(); residueNumber++) {
            String threeLetterCode = AminoAcid.Family.resolveOneLetterCode(sequence.substring(residueNumber - 1, residueNumber)).getThreeLetterCode();
            AminoAcid aminoAcid = (AminoAcid) ProteinParser.createGroup(threeLetterCode, IdentifierFactory.createResidueIdentifier(residueNumber,  ""), false, false);
            chain.addGroup(aminoAcid);
        }
        protein.addChain(chain);
        return protein;
    }

    @Override
    public String getIdentifier() {
        return identifier;
    }

    @Override
    public String getQuerySequence() {
        return querySequence;
    }

    @Override
    public Protein getQueryProtein() {
        return queryProtein;
    }

    @Override
    public Chain getQueryChain() {
        return queryChain;
    }

    @Override
    public void setHomologousEntryContainer(UniProtHomologousEntryContainer homologousEntryContainer) {
        this.homologousEntryContainer = homologousEntryContainer;
    }

    @Override
    public UniProtHomologousEntryContainer getHomologousEntryContainer() {
        return homologousEntryContainer;
    }

    @Override
    public void setHomologousPdbChains(List<Chain> homologousPdbChains) {
        this.homologousPdbChains = homologousPdbChains;
    }

    @Override
    public List<Chain> getHomologousPdbChains() {
        return homologousPdbChains;
    }

    @Override
    public Protein getReferenceProtein() {
        return referenceProtein;
    }

    @Override
    public void setReferenceProtein(Protein referenceProtein) {
        this.referenceProtein = referenceProtein;
    }

    @Override
    public Chain getReferenceChain() {
        return referenceChain;
    }

    @Override
    public void setReferenceChain(Chain referenceChain) {
        this.referenceChain = referenceChain;
    }

    @Override
    public Map<String, String> getAlignmentMap() {
        return alignmentMap;
    }

    @Override
    public void setAlignmentMap(Map<String, String> alignmentMap) {
        this.alignmentMap = alignmentMap;
    }

    @Override
    public List<SequenceConservationProfile> getSequenceConservationProfile() {
        return referenceChain.aminoAcids()
                .map(AminoAcid::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(SequenceConservationProfile.class))
                .collect(Collectors.toList());
    }

    @Override
    public boolean predictMutationEffect(int position, AminoAcid.Family targetAminoAcid) {
        AminoAcid originalAminoAcid = (AminoAcid) queryChain.getGroups().get(position - 1);
        int renumberedPosition = originalAminoAcid.getResidueIdentifier().getResidueNumber();
        logger.info("evaluating mutation effect {}->{} at position {} (renumbered: {})",
                originalAminoAcid.getOneLetterCode(),
                targetAminoAcid.getOneLetterCode(),
                position,
                renumberedPosition);

        // create mutated protein
        Protein mutatedProtein = mutateResidue(position, targetAminoAcid);
        Chain mutatedChain = mutatedProtein.select()
                .chainName(referenceChain.getChainIdentifier().getChainId())
                .asChain();
        AminoAcid mutatedAminoAcid = mutatedChain.select()
                .residueNumber(position)
                .asAminoAcid();

        // extract environments
        List<List<Group>> originalEnvironments = extractStructureFragments(homologousPdbChains, renumberedPosition);
        List<Group> originalEnvironment = extractStructureFragment(referenceChain, renumberedPosition).get();
        List<Group> mutatedEnvironment = extractStructureFragment(mutatedChain, renumberedPosition).get();

        // run mmm
//        executeMolecularMinerJob(originalEnvironments);

        FeatureVector originalEnvironmentFeatureVector = new FeatureVector(originalEnvironment);
        FeatureVector originalAminoAcidFeatureVector = new FeatureVector(originalAminoAcid);
        FeatureVector mutatedEnvironmentFeatureVector = new FeatureVector(mutatedEnvironment);
        FeatureVector mutatedAminoAcidFeatureVector = new FeatureVector(mutatedAminoAcid);
        FeatureVector environmentDelta = new FeatureVector(originalEnvironmentFeatureVector, mutatedEnvironmentFeatureVector);
        FeatureVector aminoAcidDelta = new FeatureVector(originalAminoAcidFeatureVector, mutatedAminoAcidFeatureVector);

        //TODO resume by creating full-fledged feature vectors

        return false;
    }

    private Protein mutateResidue(int renumberedPosition, AminoAcid.Family targetAminoAcid) {
        ResidueMutatorServiceImpl residueMutatorService = new ResidueMutatorServiceImpl();
        Protein mutatedProtein = residueMutatorService.mutateResidue(referenceProtein, referenceChain.getChainIdentifier().getChainId(), renumberedPosition, targetAminoAcid);
        CommonFeatureAnnotator.annotateProtein(mutatedProtein);
        return mutatedProtein;
    }

    private List<List<Group>> extractStructureFragments(List<Chain> chains, int renumberedPosition) {
        logger.info("extracting structural environments {} A around renumbered position {}", INTERACTION_CUTOFF, renumberedPosition);
        return chains.stream()
                .map(chain -> extractStructureFragment(chain, renumberedPosition))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }

    private Optional<List<Group>> extractStructureFragment(Chain chain, int renumberedPosition) {
        try {
            Atom referenceAtom = chain.select()
                    .residueNumber(renumberedPosition)
                    .asAminoAcid()
                    .getCa();
            List<Group> groups = chain.select()
                    .groupDistance(referenceAtom, INTERACTION_CUTOFF)
                    .asFilteredGroups()
                    .collect(Collectors.toList());
            return Optional.of(groups);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private void executeMolecularMinerJob(List<List<Group>> environments) {
        Path basePath = MOLECULAR_MINER_WORKING_PATH.resolve(identifier);
        Path structurePath = basePath.resolve("structures");
        Path outputPath = basePath.resolve("mmm-out");
        logger.info("writing structural fragments to {}", structurePath);

        // ensure directories are created and empty
        try {
            if(Files.exists(basePath)) {
                deleteDirectory(basePath);
            }
            Files.createDirectories(basePath);
            if (!Files.exists(structurePath)) {
                Files.createDirectories(structurePath);
            }
            if(!Files.exists(outputPath)) {
                Files.createDirectories(outputPath);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }

        // write all fragments to tmp dir
        environments.forEach(environment -> {
            try {
                ChainIdentifier chainIdentifier = environment.get(0).getParentChain().getChainIdentifier();
                String pdbId = chainIdentifier.getFullName();
                Files.write(structurePath.resolve(pdbId + ".pdb"),
                        ("HEADER    UNKNOWN PROTEIN                         01-JAN-00   " +
                                pdbId.substring(0, 4).toUpperCase() + System.lineSeparator() +
                                environment.stream()
                                        .map(Group::getPdbRepresentation)
                                        .collect(Collectors.joining())).getBytes());
            } catch (IOException e) {
                e.printStackTrace();
                throw new UncheckedIOException(e);
            }
        });
        logger.info("wrote {} environments", environments.size());

        // start mining process
        try {
            logger.info("starting molecular miner thread");
            MolecularMinerBridge.submitJob(structurePath, outputPath).get();
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    public void deleteDirectory(Path directory) throws IOException {
        Files.walkFileTree(directory, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.delete(file);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                Files.delete(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }
}
