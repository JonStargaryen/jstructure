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

    private final String identifier;
    private final String querySequence;
    private final Protein protein;
    private final Chain chain;
    private UniProtHomologousEntryContainer homologousEntryContainer;
    private List<Chain> homologousPdbChains;

    MutationEffectPredictionImpl(String identifier, String querySequence) {
        this.identifier = identifier;
        this.querySequence = querySequence;
        this.protein = createProtein(identifier, querySequence);
        this.chain = protein.getChains().get(0);
    }

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
    public Protein getProtein() {
        return protein;
    }

    @Override
    public Chain getChain() {
        return chain;
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
