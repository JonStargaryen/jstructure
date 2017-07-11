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
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtHit;

import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The implementation of the mutation effect prediction pipeline.
 * Created by bittrich on 7/11/17.
 */
public class MutationEffectPredictionServiceImpl implements MutationEffectPredictionService {
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
    public MutationEffectPrediction getMutationEffectPrediction(String identifier, String sequence) {
        MutationEffectPrediction mutationEffectPrediction = new MutationEffectPredictionImpl(identifier, sequence);

        executeBlastQuery(mutationEffectPrediction);

        return mutationEffectPrediction;
    }

    void executeBlastQuery(MutationEffectPrediction mutationEffectPrediction) {
        Protein protein = mutationEffectPrediction.getProtein();
        Chain chain = mutationEffectPrediction.getChain();

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

        // renumber everything (sequences + structures) for a consistent numbering
        List<String> sequences = Stream.concat(homologousEntryContainer.getUniProtHits()
                .stream()
                .map(UniProtHit::getEntry)
                .map(entry -> ">" + entry.getPrimaryUniProtAccession().getValue() + System.lineSeparator() + entry.getSequence().getValue()),
                homologousChains.stream()
                .map(c -> ">" + c.getChainIdentifier().getFullName() + System.lineSeparator() + c.getAminoAcidSequence()))
                .collect(Collectors.toList());
        try {
            String alignmentResult = clustalOmegaQuery.process(sequences);
            System.out.println(alignmentResult);
        } catch (ExecutionException e) {
            e.printStackTrace();
        }
    }

    private Protein createProtein(ProteinIdentifier proteinIdentifier) {
        Protein protein = ProteinParser.localPdb(proteinIdentifier.getPdbId()).minimalParsing(true).parse();

        accessibleSurfaceAreaCalculator.process(protein);
        energyProfileCalculator.process(protein);
        loopFractionCalculator.process(protein);
        plipAnnotator.process(protein);

        return protein;
    }
}
