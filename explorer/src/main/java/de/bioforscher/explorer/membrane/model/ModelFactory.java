package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.clustalo.ClustalOmegaQuery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Construct an model instance.
 * Created by bittrich on 3/17/17.
 */
public class ModelFactory {
    private static final Logger logger = LoggerFactory.getLogger(ModelFactory.class);

    public static List<ExplorerChain> createChains(String representativeChainId, List<String> chainIds) {
        // all proteins - key: pdbId, value: instance of the whole protein
        Map<String, Protein> proteins = chainIds.stream()
                .map(chainId -> chainId.split("\\.")[0])
                .limit(10)
                .distinct()
                .peek(pdbId -> logger.info("[{}] fetching homologous protein {}", representativeChainId, pdbId))
                .collect(Collectors.toMap(Function.identity(), pdbId -> ProteinParser.source(pdbId).parse()));

        return chainIds.stream()
                .limit(10)
                .map(chainId -> mapToExplorerChain(chainId, proteins, representativeChainId, chainIds))
                .collect(Collectors.toList());
    }

    public static MultiSequenceAlignment createMultiSequenceAlignment(String representativeChainId, List<ExplorerChain> explorerChains) {
        String alignment = new ClustalOmegaQuery().process(explorerChains.stream()
                .map(ExplorerChain::getSequence)
                .collect(Collectors.toList()));
        return new MultiSequenceAlignment(representativeChainId, alignment);
    }

    private static ExplorerChain mapToExplorerChain(String chainId, Map<String, Protein> proteins, String representativeChainId, List<String> homologous) {
        String[] split = chainId.split("\\.");
        Protein protein = proteins.get(split[0]);
        return new ExplorerChain(protein.select().chainName(split[1]).asChain(), homologous, representativeChainId);
    }
}
