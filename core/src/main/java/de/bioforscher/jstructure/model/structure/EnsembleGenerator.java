package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.identifier.PdbChainId;
import de.bioforscher.jstructure.model.identifier.PdbId;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Creates collections of proteins or protein chains based on some unifying characteristic. Support are e.g. sequence
 * clusters from the PDB or UniProt as well as Pfam data.
 * Created by bittrich on 4/27/17.
 */
public class EnsembleGenerator {
    private static final Logger logger = LoggerFactory.getLogger(EnsembleGenerator.class);

    //TODO some builder patterns to create features on-the-fly

    public enum PdbSequenceClusterCutoff {
        CUTOFF100(100),
        CUTOFF95(95),
        CUTOFF90(90),
        CUTOFF70(70),
        CUTOFF50(50),
        CUTOFF40(40),
        CUTOFF30(30);

        private int cutoff;

        PdbSequenceClusterCutoff(int cutoff) {
            this.cutoff = cutoff;
        }

        public int getCutoff() {
            return cutoff;
        }
    }

    public static Stream<Protein> createSequenceClusterFromPdb(PdbChainId referenceId, PdbSequenceClusterCutoff cutoff) {
        return createStreamFromPdbChainIds(referenceId.getFullName(), PDBEnsembleQuery.fetchSequenceCluster(referenceId, cutoff));
    }

    public static Stream<Protein> createSequenceClusterFromUniProt() {
        //TODO impl
        throw new UnsupportedOperationException("not implemented!");
    }

    public static Stream<Protein> createSequenceClusterFromPfam(String pfamId) {
        return createStreamFromPdbChainIds(pfamId, new SiftsMappingProvider().mapPfamIdToPdbIds(pfamId));
    }

    private static Stream<Protein> createStreamFromPdbChainIds(String referenceId, List<PdbChainId> pdbChainIds) {
        logger.info("[] cluster consists of []", referenceId, pdbChainIds);
        return pdbChainIds.parallelStream()
                .map(pdbChainId -> ProteinParser.source(pdbChainId.getPdbId().getPdbId()).parse().select().chainName(pdbChainId.getChainId()).asChainContainer())
                .map(Protein.class::cast);
    }

    private static Protein mapToProtein(PdbChainId pdbChainId) {
        return (Protein) ProteinParser.source(pdbChainId.getPdbId().getPdbId()).parse().select().chainName(pdbChainId.getChainId()).asChainContainer();
    }


    static class PDBEnsembleQuery {
        static final String CLUSTER_FETCH_URL = "http://www.rcsb.org/pdb/rest/sequenceCluster?cluster=%d&structureId=%s.%s";

        /**
         * Queries the PDB for all chains which share a given sequence similarity to a reference chain.
         * @param pdbChainId the reference pdb id
         * @param cutoff the desired sequence similarity threshold
         * @return all chain ids (in the format pdbId_chainId) which are similar to the reference - they are sorted, the
         *      first entry is the most relevant one
         */
        static List<PdbChainId> fetchSequenceCluster(PdbChainId pdbChainId, PdbSequenceClusterCutoff cutoff) {
            try {
                logger.debug("fetching {}-cluster fom PDB for {}_{}", cutoff.getCutoff(), pdbChainId);
                Document document = Jsoup.connect(String.format(CLUSTER_FETCH_URL, cutoff.getCutoff(), pdbChainId.getPdbId().getPdbId(), pdbChainId.getChainId())).get();
                return document.getElementsByTag("pdbChain").stream()
                        .map(pdbChain -> pdbChain.attr("name"))
                        .map(pdbChain -> pdbChain.split("\\."))
                        .map(split -> PdbChainId.createFromChainId(PdbId.createFromPdbId(split[0]), split[1]))
                        .collect(Collectors.toList());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }
}
