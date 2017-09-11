package de.bioforscher.jstructure.membrane;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Query the PDBTM and obtain the current list of non-redundant beta-barrel chain ids. Ensure that the corresponding
 * chain is present in the PDB, otherwise the chain is ignored and will not be present in the local list.
 */
public class T02_DownloadPdbtmBetaNrDataset {
    private static final Logger logger = LoggerFactory.getLogger(T02_DownloadPdbtmBetaNrDataset.class);
    private static final String FETCH_URL = "http://pdbtm.enzim.hu/data/pdbtm_beta_nr.list";

    public static void main(String[] args) {
        String selectedIds = MembraneConstants.lines(FETCH_URL)
                .filter(T02_DownloadPdbtmBetaNrDataset::handlePdbId)
                .collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(MembraneConstants.PDBTM_BETA_DATASET_DIRECTORY.resolve("pdbtm_beta_nr.list"), selectedIds);
    }

    private static boolean handlePdbId(String chainId) {
        try {
            String pdbId = chainId.split("_")[0];
            String chainIdentifier = chainId.split("_")[1];
            Structure structure = StructureParser.source(pdbId).minimalParsing(true).parse();
            Optional<Chain> chain = structure.select()
                    .chainId(chainIdentifier)
                    .asOptionalChain();
            if(!chain.isPresent()) {
                logger.warn("{} is missing in structure", chainId);
                return false;
            }
            MembraneConstants.write(MembraneConstants.PDBTM_BETA_DATASET_PDB_DIRECTORY.resolve(pdbId + ".pdb"),
                    structure.getPdbRepresentation());
            logger.info("{} is present", chainId);
            return true;
        } catch (Exception e) {
            logger.warn("{} is obsolete", chainId);
            return false;
        }
    }
}
