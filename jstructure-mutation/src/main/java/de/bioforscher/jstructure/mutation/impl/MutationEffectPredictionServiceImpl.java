package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.align.ClustalOmegaRestQuery;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologyAnnotator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.mutation.MutationEffectPrediction;
import de.bioforscher.jstructure.mutation.MutationEffectPredictionService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtHit;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * The implementation of the mutation effect prediction pipeline.
 * Created by bittrich on 7/11/17.
 */
public class MutationEffectPredictionServiceImpl implements MutationEffectPredictionService {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionServiceImpl.class);
    private final UniProtHomologyAnnotator uniProtHomologyAnnotator;
    private final AccessibleSurfaceAreaCalculator accessibleSurfaceAreaCalculator;
    private final EnergyProfileCalculator energyProfileCalculator;
    private final PLIPAnnotator plipAnnotator;
    private final LoopFractionCalculator loopFractionCalculator;
    private final ClustalOmegaRestQuery clustalOmegaQuery;

    public MutationEffectPredictionServiceImpl() {
        this.uniProtHomologyAnnotator = new UniProtHomologyAnnotator();
        this.accessibleSurfaceAreaCalculator = new AccessibleSurfaceAreaCalculator();
        this.loopFractionCalculator = new LoopFractionCalculator();
        this.plipAnnotator = new PLIPAnnotator();
        this.energyProfileCalculator = new EnergyProfileCalculator();
        this.clustalOmegaQuery = new ClustalOmegaRestQuery();

        //TODO unify behaviour - global config?
        ProteinParser.OptionalSteps.setLocalPdbDirectory(Paths.get("/home/bittrich/pdb/"));
    }

    @Override
    public MutationEffectPrediction getMutationEffectPrediction(String identifier, String sequence) throws ExecutionException {
        logger.info("starting new job '{}' with sequence {}", identifier, sequence);
        MutationEffectPrediction mutationEffectPrediction = new MutationEffectPredictionImpl(identifier, sequence);

        createMultiSequenceAlignment(mutationEffectPrediction);

        return mutationEffectPrediction;
    }

    /**
     * Create a multi-sequence alignment with the given query sequence. Will run BLAST against SWISS-PROT and retrieve
     * all UniProt sequences returned. Of those entries all associated protein chains in the PDB will be fetched and
     * annotated with features. All those sequence (query, UniProt sequences and those of the PDB chains) will lastly be
     * aligned using ClustalOmega to achieve a unified numbering of all residues.
     * @param mutationEffectPrediction the container to manipulate
     */
    void createMultiSequenceAlignment(MutationEffectPrediction mutationEffectPrediction) throws ExecutionException {
        // create wrapping pseudo-instances
        Protein protein = mutationEffectPrediction.getQueryProtein();
        Chain chain = mutationEffectPrediction.getQueryChain();

        // execute BLAST
        uniProtHomologyAnnotator.process(protein);

        // link associated UniProt entries
        UniProtHomologousEntryContainer homologousEntryContainer = chain.getFeatureContainer().getFeature(UniProtHomologousEntryContainer.class);
        mutationEffectPrediction.setHomologousEntryContainer(homologousEntryContainer);

        // retrieve homologous protein structures
        Map<ProteinIdentifier, Protein> proteinMap = homologousEntryContainer.getHomologousChains()
                .stream()
                .map(ChainIdentifier::getProteinIdentifier)
                .distinct()
                .collect(Collectors.toMap(Function.identity(), this::createProtein));
        List<Chain> homologousChains = homologousEntryContainer.getHomologousChains()
                .stream()
                .map(chainIdentifier -> proteinMap.get(chainIdentifier.getProteinIdentifier()).select().chainName(chainIdentifier.getChainId()).asChain())
                .collect(Collectors.toList());
        mutationEffectPrediction.setHomologousPdbChains(homologousChains);

        // execute ClustalOmega on all sequences
        List<String> sequences = new ArrayList<>();
        sequences.add(">query" + System.lineSeparator() + mutationEffectPrediction.getQuerySequence());
        homologousEntryContainer.getUniProtHits()
                .stream()
                .map(UniProtHit::getEntry)
                .map(entry -> ">" + entry.getPrimaryUniProtAccession().getValue() + System.lineSeparator() + entry.getSequence().getValue())
                .forEach(sequences::add);
        homologousChains.stream()
                .map(c -> ">" + c.getChainIdentifier().getFullName() + System.lineSeparator() + c.getAminoAcidSequence())
                .forEach(sequences::add);
        String alignmentResult = clustalOmegaQuery.process(sequences);

        // renumber everything with respect to the query sequence
        //TODO impl
    }

    /**
     * Create a Protein instance and compute all necessary features.
     * @param proteinIdentifier the identifier to process
     * @return the created and annotated instance
     */
    private Protein createProtein(ProteinIdentifier proteinIdentifier) {
        Protein protein = ProteinParser.localPdb(proteinIdentifier.getPdbId()).minimalParsing(true).parse();

        accessibleSurfaceAreaCalculator.process(protein);
        energyProfileCalculator.process(protein);
        loopFractionCalculator.process(protein);
        plipAnnotator.process(protein);

        return protein;
    }
}
