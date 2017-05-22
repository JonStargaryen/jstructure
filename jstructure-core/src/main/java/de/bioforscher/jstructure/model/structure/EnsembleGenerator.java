package de.bioforscher.jstructure.model.structure;

/**
 * Creates collections of proteins or protein chains based on some unifying characteristic. Support are e.g. sequence
 * clusters from the PDB or UniProt as well as Pfam data.
 * Created by bittrich on 4/27/17.
 */
public class EnsembleGenerator {
    //TODO reimpl
//    private static final Logger logger = LoggerFactory.getLogger(EnsembleGenerator.class);
//    private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();
//
//    public enum PdbSequenceClusterCutoff {
//        CUTOFF100(100),
//        CUTOFF95(95),
//        CUTOFF90(90),
//        CUTOFF70(70),
//        CUTOFF50(50),
//        CUTOFF40(40),
//        CUTOFF30(30);
//
//        private int cutoff;
//
//        PdbSequenceClusterCutoff(int cutoff) {
//            this.cutoff = cutoff;
//        }
//
//        public int getCutoff() {
//            return cutoff;
//        }
//    }
//
//    public static List<Chain> createChainsFromPdbSequenceCluster(ChainIdentifier referenceId, PdbSequenceClusterCutoff cutoff) {
//        return createStreamFromPdbChainIds(referenceId.getFullName(), PDBEnsembleQuery.fetchSequenceCluster(referenceId, cutoff));
//    }
//
//    public static List<Chain> createChainsFromUniProtSequenceCluster() {
//        //TODO impl
//        throw new UnsupportedOperationException("not implemented!");
//    }
//
//    public static List<Chain> createChainsFromPfamCluster(String pfamId) {
//        return createStreamFromPdbChainIds(pfamId, new SiftsMappingProvider().mapPfamIdToPdbIds(pfamId));
//    }
//
//    private static List<Chain> createStreamFromPdbChainIds(String referenceId, List<ChainIdentifier> pdbChainIds) {
//        logger.info("[{}] cluster consists of {}", referenceId, pdbChainIds);
//
//        Map<ProteinIdentifier, Protein> uniqueProteins = pdbChainIds.parallelStream()
//                .map(ChainIdentifier::getPdbId)
//                .distinct()
//                .peek(pdbId -> logger.info("[{}] fetching pdb structure {}", referenceId, pdbId))
//                .collect(Collectors.toMap(Function.identity(), pdbId -> ProteinParser.source(pdbId.getPdbId()).parse()));
//
//        List<Chain> chains = pdbChainIds.stream()
//                .map(pdbChainId -> uniqueProteins.get(pdbChainId.getPdbId()).select().chainName(pdbChainId.getChainId()).asChain())
//                .collect(Collectors.toList());
//
//        chains.parallelStream().forEach(plipAnnotator::process);
//
//        // TODO builder pattern for features, alignments etc
//        StructureAligner.align(chains);
////        SequenceAligner.align(chains);
//
//        handleUniProtInformation(chains, referenceId);
//
//        return chains;
//    }
//
//    private static void handleUniProtInformation(List<Chain> chains, String referenceId) {
//        //TODO reimpl
////        chains.forEach(new SiftsMappingProvider()::processUniProtId);
////        List<String> uniProtIds = chains.stream()
////                .map(chain -> chain.getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING))
////                .distinct()
////                .map(ChainSiftsMapping::getUniProtId)
////                .distinct()
////                .collect(Collectors.toList());
////
////        logger.info("[{}] associated UniProt ids: {}", referenceId, uniProtIds);
////
////        Map<String, UniProtAnnotationContainer> uniProtAnnotationContainers = uniProtIds.stream()
////                .peek(uniProtId -> logger.info("[{}] fetching UniProt information for {}", referenceId, uniProtId))
////                .collect(Collectors.toMap(Function.identity(), uniProtId -> new UniProtAnnotator().processUniProtId(uniProtId)));
////
////        chains.forEach(chain -> {
////            String uniProtId = chain.getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING).getUniProtId();
////            chain.setFeature(UniProtAnnotator.UNIPROT_ANNOTATION, uniProtAnnotationContainers.get(uniProtId));
////        });
//    }
//
//    static class StructureAligner {
//        private static final SVDSuperimposer svdSuperimposer = new SVDSuperimposer();
//        private Chain reference;
//        private List<Chain> chains;
//
//        private StructureAligner(List<Chain> chains) {
//            this.reference = chains.get(0);
//            this.chains = chains.subList(1, chains.size());
//        }
//
//        static void align(List<Chain> chains) {
//            StructureAligner structureAligner = new StructureAligner(chains);
//            structureAligner.alignInternally();
//        }
//
//        private void alignInternally() {
//            chains.forEach(this::alignInternally);
//        }
//
//        private void alignInternally(Chain chain) {
//            StructureAlignmentResult alignmentResult = svdSuperimposer.align(reference, chain);
//            alignmentResult.transform(chain);
//            //TODO reimpl
////            chain.setFeature(SVDSuperimposer.TRANSFORMATION, new LinearAlgebraAtom.Transformation(alignmentResult.getTranslation(), alignmentResult.getRotation()));
//        }
//    }
//
//    static class SequenceAligner {
//        //TODO impl
//    }
//
//    static class PDBEnsembleQuery {
//        static final String CLUSTER_FETCH_URL = "http://www.rcsb.org/pdb/rest/sequenceCluster?cluster=%d&structureId=%s.%s";
//
//        /**
//         * Queries the PDB for all chains which share a given sequence similarity to a reference chain.
//         * @param pdbChainId the reference pdb id
//         * @param cutoff the desired sequence similarity threshold
//         * @return all chain ids (in the format pdbId_chainId) which are similar to the reference - they are sorted, the
//         *      first entry is the most relevant one
//         */
//        static List<ChainIdentifier> fetchSequenceCluster(ChainIdentifier pdbChainId, PdbSequenceClusterCutoff cutoff) {
//            try {
//                logger.debug("fetching {}-cluster fom PDB for {}", cutoff.getCutoff(), pdbChainId);
//                Document document = Jsoup.connect(String.format(CLUSTER_FETCH_URL, cutoff.getCutoff(), pdbChainId.getPdbId().getPdbId(), pdbChainId.getChainId())).get();
//                return document.getElementsByTag("pdbChain").stream()
//                        .map(pdbChain -> pdbChain.attr("name"))
//                        .map(pdbChain -> pdbChain.split("\\."))
//                        .map(split -> ChainIdentifier.createFromChainId(ProteinIdentifier.createFromPdbId(split[0]), split[1]))
//                        .collect(Collectors.toList());
//            } catch (IOException e) {
//                throw new UncheckedIOException(e);
//            }
//        }
//    }
}
