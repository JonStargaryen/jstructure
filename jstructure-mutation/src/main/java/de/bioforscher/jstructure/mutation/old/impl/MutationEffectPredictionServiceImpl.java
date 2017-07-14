package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.align.impl.ClustalOmegaRestQuery;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologyAnnotator;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.mutation.old.MutationEffectPredictionService;
import de.bioforscher.jstructure.mutation.old.MutationJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtHit;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

/**
 * The implementation of the mutation effect prediction pipeline.
 * Created by bittrich on 7/11/17.
 */
@Deprecated
public class MutationEffectPredictionServiceImpl implements MutationEffectPredictionService {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionServiceImpl.class);
    private static final String QUERY_ID = "query";
    private final UniProtHomologyAnnotator uniProtHomologyAnnotator;
    private final ClustalOmegaRestQuery clustalOmegaQuery;

    public MutationEffectPredictionServiceImpl() {
        this.uniProtHomologyAnnotator = new UniProtHomologyAnnotator();
        this.clustalOmegaQuery = new ClustalOmegaRestQuery();

        //TODO unify behaviour - global config?
        ProteinParser.OptionalSteps.setLocalPdbDirectory(Paths.get("/home/bittrich/pdb/"));
    }

    @Override
    public MutationJob createMutationJob(String identifier, String sequence) throws ExecutionException {
        logger.info("starting new job '{}' with sequence {}", identifier, sequence);
        MutationJob mutationJob = new MutationJobImpl(identifier, sequence);

        createMultiSequenceAlignment(mutationJob);

        return mutationJob;
    }

    /**
     * Create a multi-sequence alignment with the given query sequence. Will run BLAST against SWISS-PROT and retrieve
     * all UniProt sequences returned. Of those entries all associated protein chains in the PDB will be fetched and
     * annotated with features. All those sequence (query, UniProt sequences and those of the PDB chains) will lastly be
     * aligned using ClustalOmega to achieve a unified numbering of all residues.
     * @param mutationJob the container to manipulate
     */
    void createMultiSequenceAlignment(MutationJob mutationJob) {
        // create wrapping pseudo-instances
        Protein queryProtein = mutationJob.getQueryProtein();
        Chain queryChain = mutationJob.getQueryChain();

        // execute BLAST
        uniProtHomologyAnnotator.process(queryProtein);

        // link associated UniProt entries
        UniProtHomologousEntryContainer homologousEntryContainer = queryChain.getFeatureContainer().getFeature(UniProtHomologousEntryContainer.class);
        mutationJob.setHomologousEntryContainer(homologousEntryContainer);

        // retrieve homologous protein structures
        Map<ProteinIdentifier, Protein> proteinMap = homologousEntryContainer.getHomologousChains()
                .stream()
                .map(ChainIdentifier::getProteinIdentifier)
                .distinct()
                .map(proteinIdentifier -> new Pair<>(proteinIdentifier, createProtein(proteinIdentifier)))
                .filter(pair -> pair.getRight().isPresent())
                .map(pair -> new Pair<>(pair.getLeft(), pair.getRight().get()))
                .collect(Collectors.toMap(Pair::getLeft, Pair::getRight));
        // if structure was provided, annotate it with features
        if(mutationJob.getReferenceChain() != null) {
            CommonFeatureAnnotator.annotateProtein(mutationJob.getReferenceProtein());
        }
        List<Chain> homologousChains = homologousEntryContainer.getHomologousChains()
                .stream()
                .map(chainIdentifier -> proteinMap.get(chainIdentifier.getProteinIdentifier()).select().chainName(chainIdentifier.getChainId()).asChain())
                .collect(Collectors.toList());
        mutationJob.setHomologousPdbChains(homologousChains);
        logger.info("{} structure hits: {}", homologousChains.size(), homologousChains.stream().map(Chain::getChainIdentifier).collect(Collectors.toList()));

        // execute ClustalOmega on all sequences
        List<String> sequences = new ArrayList<>();
        sequences.add(">" + QUERY_ID + System.lineSeparator() + mutationJob.getQuerySequence());
        homologousEntryContainer.getUniProtHits()
                .stream()
                .map(UniProtHit::getEntry)
                .map(entry -> ">" + entry.getPrimaryUniProtAccession().getValue() + System.lineSeparator() + entry.getSequence().getValue())
                .forEach(sequences::add);
        homologousChains.stream()
                .map(chain -> ">" + chain.getChainIdentifier().getFullName() + System.lineSeparator() + chain.getAminoAcidSequence())
                .forEach(sequences::add);
        Map<String, String> alignmentMap = clustalOmegaQuery.align(sequences).getAlignedSequences();

        // renumber everything with respect to the query sequence
        renumber(queryChain, alignmentMap.get(QUERY_ID));
        // if structure was provided, renumber that too
        if(mutationJob.getReferenceChain() != null) {
            renumber(mutationJob.getReferenceChain(), alignmentMap.get(QUERY_ID));
        }
        mutationJob.getHomologousPdbChains().forEach(chain -> renumber(chain, alignmentMap.get(chain.getChainIdentifier().getFullName())));

        // create sequence conservation profile for each residue
        SequenceConservationAnnotator.process(queryChain, alignmentMap);

        // assign reference structure
        if(homologousChains.size() > 0) {
            // only needed when based on sequence and no structure is available
            if(mutationJob.getReferenceChain() == null) {
                // TODO select more reasonable
                Chain referenceChain = homologousChains.get(0);
                mutationJob.setReferenceChain(referenceChain);
                mutationJob.setReferenceProtein(referenceChain.getParentProtein());
            }
        } else {
            logger.warn("no structure hits for {}", mutationJob.getIdentifier());
        }
    }

    private void renumber(Chain chain, String alignmentString) {
        int alignmentLength = alignmentString.length();
        int consumedAminoAcidsInChain = 0;
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        for(int sequencePosition = 0; sequencePosition < alignmentLength; sequencePosition++) {
            char characterInAlignment = alignmentString.charAt(sequencePosition);
            if(characterInAlignment == '-') {
                continue;
            }
            AminoAcid aminoAcid = aminoAcids.get(consumedAminoAcidsInChain);
            consumedAminoAcidsInChain++;
            // renumber to position in alignment string (+1 for classic Java offset)
            aminoAcid.setResidueIdentifier(IdentifierFactory.createResidueIdentifier(sequencePosition + 1));
        }
    }

    /**
     * Create a Protein instance and compute all necessary features.
     * @param proteinIdentifier the identifier to process
     * @return the created and annotated instance
     */
    private Optional<Protein> createProtein(ProteinIdentifier proteinIdentifier) {
        try {
            Protein protein = ProteinParser.source(proteinIdentifier.getPdbId())
                    .minimalParsing(true)
                    .parse();
            CommonFeatureAnnotator.annotateProtein(protein);
            return Optional.of(protein);
        } catch (Exception e) {
            return Optional.empty();
        }
    }
}
