package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.clustalo.ClustalOmegaRestQuery;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Construct an model instance.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerModelFactory {
    private static final Logger logger = LoggerFactory.getLogger(ExplorerModelFactory.class);
    //TODO some 'worker'-patterns
    private static final AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolve(PLIPAnnotator.PLIP_INTERACTIONS);
    private static final AbstractFeatureProvider uniprotAnnotator = FeatureProviderRegistry.resolve(UniProtAnnotator.UNIPROT_ANNOTATION);
    static final String TRANSFORMATION = "TRANSFORMATION";

    public static ExplorerCluster createCluster(String representativeChainId,
                                                  List<String> chainIds) {
        // all proteins - key: pdbId, value: instance of the whole protein
        Map<String, Protein> proteins = chainIds.stream()
                .map(chainId -> chainId.split("_")[0].toLowerCase())
                .distinct()
                .peek(pdbId -> logger.info("[{}] fetching homologous protein {}", representativeChainId, pdbId))
                .collect(Collectors.toMap(Function.identity(),
                        pdbId -> ProteinParser.source(pdbId).forceProteinName(pdbId).parse()));

        proteins.values().parallelStream().forEach(uniprotAnnotator::process);
//        proteins.values().parallelStream().forEach(plipAnnotator::process);

        // the reference of the pv_structure alignment
        Chain referenceChain = proteins.get(representativeChainId.split("_")[0])
                .select()
                .chainName(representativeChainId.split("_")[1])
                .nameContainer(representativeChainId)
                .asChain();

        // extract wanted original chains
        List<Chain> chains = chainIds.stream()
                .map(chainId -> proteins.get(chainId.split("_")[0])
                        .select()
                        .chainName(chainId.split("_")[1])
                        .nameContainer(chainId)
                        .asChain())
                .collect(Collectors.toList());

        chains.forEach(chain -> {
            StructureAlignmentResult alignmentResult = new SVDSuperimposer().align(referenceChain, chain);
            logger.info("[{}] rmsd of {} to reference is {}", representativeChainId, getGlobalId(chain), alignmentResult.getAlignmentScore());
            alignmentResult.transform(chain);
            chain.setFeature(TRANSFORMATION, new LinearAlgebraAtom.Transformation(alignmentResult.getTranslation(), alignmentResult.getRotation()));
        });

        logger.info("[{}] creating multi-sequence alignment by clustal omega", representativeChainId);
        //TODO case of no homologous
        // align sequences
        List<String> sequences = chains.stream()
                .map(chain -> ">" + getGlobalId(chain) + "|PDBID|CHAIN|SEQUENCE" + System.lineSeparator() + chain.getAminoAcidSequence())
                .collect(Collectors.toList());
        try {
            String alignmentString = new ClustalOmegaRestQuery().process(sequences);

            // split into strings describing on entry of the alignment each
            List<String> splitAlignmentString = Pattern.compile(">")
                    .splitAsStream(alignmentString)
                    .skip(1)
                    .collect(Collectors.toList());

            // key: id, value: alignedSequence with gaps
            Map<String, String> alignedSequenceMap = splitAlignmentString.stream()
                    // split lines into pairs
                    .map(alignmentLine -> alignmentLine.split("SEQUENCE"))
                    // substring for id, remove line-breaks
                    .collect(Collectors.toMap(split -> split[0].substring(0, 6),
                            split -> split[1].replaceAll("\\s", "")));

            // renumber original chains according to alignment
            chains.forEach(chain -> {
                String alignedSequence = alignedSequenceMap.get(getGlobalId(chain));
                Iterator<Group> aminoAcids = chain.select().aminoAcids().asFilteredGroups().iterator();
                for(int pos = 0; pos < alignedSequence.length(); pos++) {
                    if(alignedSequence.charAt(pos) != '-')  {
                        aminoAcids.next().setResidueNumber(pos);
                    }
                }
            });

            // create renumbered explorer chain objects which will provide additional, chain-specific information and the pdb record
            List<ExplorerChain> explorerChains = chains.stream()
                    .map(chain -> new ExplorerChain(chain, representativeChainId))
                    .collect(Collectors.toList());

            ExplorerAlignment alignment = new ExplorerAlignment(representativeChainId, chains, alignedSequenceMap);

            return new ExplorerCluster(alignment, explorerChains);
        } catch (ExecutionException e) {
            //TODO error-handling
            throw new RuntimeException(e);
        }
    }

    public static String getGlobalId(Chain chain) {
        return chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId();
    }
}
