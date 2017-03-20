package de.bioforscher.explorer.membrane.model;

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

        proteins.values().parallelStream().forEach(plipAnnotator::process);

        return chainIds.stream()
                .map(chainId -> mapToExplorerChain(chainId, proteins, representativeChainId, chainIds))
                .collect(Collectors.toList());
    }

    public static MultiSequenceAlignment createMultiSequenceAlignment(String representativeChainId,
                                                                      List<ExplorerChain> explorerChains) {
        try {
            String alignment = new ClustalOmegaRestQuery().process(explorerChains.stream()
                    .map(explorerChain -> ">" + explorerChain.getId() + "|PDBID|CHAIN|SEQUENCE" + System.lineSeparator() + explorerChain.getSequence())
                    .collect(Collectors.toList()));
            return new MultiSequenceAlignment(representativeChainId, alignment);
        } catch (ExecutionException e) {
            //TODO error-handling
            throw new RuntimeException(e);
        }
    }

    private static ExplorerChain mapToExplorerChain(String chainId,
                                                    Map<String, Protein> proteins,
                                                    String representativeChainId,
                                                    List<String> homologous) {
        String[] split = chainId.split("_");
        Protein protein = proteins.get(split[0]);
        Chain chain = protein.select().chainName(split[1]).nameContainer(chainId).asChain();

        //TODO plip interactions assignment

        return new ExplorerChain(chain, homologous, representativeChainId);
    }
}
