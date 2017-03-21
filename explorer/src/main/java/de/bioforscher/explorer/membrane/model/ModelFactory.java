package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.clustalo.ClustalOmegaRestQuery;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Construct an model instance.
 * Created by bittrich on 3/17/17.
 */
public class ModelFactory {
    private static final Logger logger = LoggerFactory.getLogger(ModelFactory.class);
    private static final AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolve(PLIPAnnotator.PLIP_INTERACTIONS);

    public static List<ExplorerChain> createChains(String representativeChainId,
                                                   List<String> chainIds) {
        // all proteins - key: pdbId, value: instance of the whole protein
        Map<String, Protein> proteins = chainIds.stream()
                .map(chainId -> chainId.split("_")[0].toLowerCase())
                .distinct()
                .peek(pdbId -> logger.info("[{}] fetching homologous protein {}", representativeChainId, pdbId))
                .collect(Collectors.toMap(Function.identity(),
                        pdbId -> ProteinParser.source(pdbId).forceProteinName(pdbId).parse()));

//        proteins.values().parallelStream().forEach(plipAnnotator::process);

        //TODO plip interactions assignment and rotating according to alignment

        // the reference of the structure alignment
        Chain referenceChain = proteins.get(representativeChainId.split("_")[0])
                .select()
                .chainName(representativeChainId.split("_")[1])
                .nameContainer(representativeChainId)
                .asChain();

        // extract wanted chains
        List<Chain> chains = chainIds.stream()
                .map(chainId -> proteins.get(chainId.split("_")[0])
                        .select()
                        .chainName(chainId.split("_")[1])
                        .nameContainer(chainId)
                        .asChain())
                .collect(Collectors.toList());

        // align chains to reference
        chains.forEach(chain -> {
            StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(referenceChain, chain);
            logger.info("[{}] rmsd of {} to reference is {}", representativeChainId, chain.getParentProtein().getName() + "_" + chain.getChainId(), alignmentResult.getAlignmentScore());
            alignmentResult.transform(chain);
        });

        return chains.stream()
                .map(chain -> new ExplorerChain(chain, chainIds, representativeChainId))
                .collect(Collectors.toList());
    }

    public static MultiSequenceAlignment createMultiSequenceAlignment(String representativeChainId,
                                                                      List<ExplorerChain> explorerChains) {
        try {
            String alignment = new ClustalOmegaRestQuery().process(explorerChains.stream()
                    .map(explorerChain -> ">" + explorerChain.getId() + "|PDBID|CHAIN|SEQUENCE" + System.lineSeparator() + explorerChain.getSequence())
                    .collect(Collectors.toList()));
            return new MultiSequenceAlignment(representativeChainId, explorerChains, alignment);
        } catch (ExecutionException e) {
            //TODO error-handling
            throw new RuntimeException(e);
        }
    }
}
