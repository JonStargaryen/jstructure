package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.identifier.PdbChainId;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.EnsembleGenerator;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Creates model instances.
 * Created by bittrich on 4/6/17.
 */
public class ModelFactory {
    private static final Logger logger = LoggerFactory.getLogger(ModelFactory.class);
   private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();

    public static List<ExplorerChain> createCluster(PdbChainId entryId) {
        List<Chain> chains = EnsembleGenerator.createChainsFromPfamCluster(new SiftsMappingProvider().mapPdbChainIdToPfam(entryId));

        chains.parallelStream().forEach(plipAnnotator::process);

        List<ExplorerChain> explorerChains = chains.parallelStream()
                .map(chain -> new ExplorerChain(chain, entryId.getFullName()))
                .collect(Collectors.toList());

        return explorerChains;
    }

//    /**
//     * Gathers all information on one UniProt entry.
//     * @param entryId a PDB-chain-id in the form pdbId_chainId
//     * @return the corresponding UniProt sequence cluster
//     */
//    public static ExplorerCluster createCluster(String entryId) {
//        // fetch the 95-cluster from the PBD
//        List<String> homologousChainIds = PDBDatabaseQuery.fetchSequenceCluster(entryId.split("_")[0], entryId.split("_")[1]);
//
//        // choose first entry as representative, sorted by PDB-criterion
//        String representativeId = homologousChainIds.get(0);
//        String representativePdbId = representativeId.split("_")[0];
//        String representativeChainId = representativeId.split("_")[1];
//
//        // load all proteins - key: pdbId, value: instance of the whole protein
//        Map<String, Protein> proteins = homologousChainIds.parallelStream()
//                .map(chainId -> chainId.split("_")[0].toLowerCase())
//                .distinct()
//                .peek(pdbId -> logger.info("[{}] fetching homologous protein {}", representativeId, pdbId))
//                .collect(Collectors.toMap(Function.identity(),
//                        pdbId -> ProteinParser.source(pdbId).forceProteinName(pdbId).parse()));
//
//        // annotate PLIP data
//        proteins.values().parallelStream().forEach(asaCalculator::process);
//        proteins.values().parallelStream().forEach(plipAnnotator::process);
//        proteins.values().parallelStream().forEach(uniProtAnnotator::process);
//
//        // map residues to UniProt indices - important: this has to happen after annotating any external data such as PLIP
//        proteins.values().parallelStream().forEach(protein -> {
//            // renumber according to UniProt
//            protein.aminoAcids().forEach(group -> {
//                ResidueSiftsMapping siftsMapping = group.getFeature(ResidueSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING);
//                group.setResidueNumber(siftsMapping.getUniProtResidueNumber());
//            });
//        });
//
//        // select reference chain for the structure alignment
//        Chain referenceChain = proteins.get(representativePdbId)
//                .select()
//                .chainName(representativeChainId)
//                .asChain();
//
//        // extract wanted original chains
//        List<Chain> chains = homologousChainIds.stream()
//                .map(chainId -> proteins.get(chainId.split("_")[0])
//                        .select()
//                        .chainName(chainId.split("_")[1])
//                        .nameContainer(chainId)
//                        .asChain())
//                .collect(Collectors.toList());
//
//        // align and assign transformation object to each chain
//        chains.forEach(chain -> {
//            StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(referenceChain, chain);
//            alignmentResult.transform(chain);
//            chain.setFeature(TRANSFORMATION, new LinearAlgebraAtom.Transformation(alignmentResult.getTranslation(), alignmentResult.getRotation()));
//        });
//
//        List<ExplorerChain> explorerChains = chains.parallelStream()
//                .map(chain -> new ExplorerChain(chain, representativeId))
//                .collect(Collectors.toList());
//
//        ExplorerAlignment alignment = new ExplorerAlignment(chains, representativeId);
//
//        return new ExplorerCluster(explorerChains, alignment);
//    }
}
