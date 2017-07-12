package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.mutation.MutationEffectPrediction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

/**
 * The implementation of the data structure.
 * Created by bittrich on 7/11/17.
 */
class MutationEffectPredictionImpl implements MutationEffectPrediction {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionImpl.class);
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

    MutationEffectPredictionImpl(String identifier, String querySequence) {
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
    private Protein createProtein(String identifier, String sequence) {
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
    public boolean predictMutationEffect(int position, AminoAcid.Family targetAminoAcid) {
        List<List<Group>> environments = extractStructureFragments(position);
        executeMolecularMinerJob(environments);
        return false;
    }

    private List<List<Group>> extractStructureFragments(int position) {
        logger.info("extracting structural environments at {} with {} A cutoff", position, INTERACTION_CUTOFF);
        return homologousPdbChains.stream()
                .map(chain -> extractStructureFragment(chain, /*position*/50))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }

    private Optional<List<Group>> extractStructureFragment(Chain chain, int position) {
        try {
            Atom referenceAtom = chain.select()
                    .residueNumber(position)
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
        // write all fragments to tmp dir
        environments.forEach(environment -> {
            try {
                if(!Files.exists(structurePath)) {
                    Files.createDirectories(structurePath);
                }
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
}
