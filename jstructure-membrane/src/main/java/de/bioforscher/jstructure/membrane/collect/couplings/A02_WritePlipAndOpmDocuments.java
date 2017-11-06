package de.bioforscher.jstructure.membrane.collect.couplings;

import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.UncheckedIOException;
import java.nio.file.Path;
import java.util.Optional;

public class A02_WritePlipAndOpmDocuments {
    private static final Logger logger = LoggerFactory.getLogger(A02_WritePlipAndOpmDocuments.class);
    private static final Path directory = MembraneConstants.COUPLING_DIRECTORY;

    public static void main(String[] args) {
        MembraneConstants.lines(directory.resolve("ids.list"))
                .forEach(A02_WritePlipAndOpmDocuments::handleLine);
    }

    private static void handleLine(String chainId) {
        try {
            String pdbId = chainId.split("_")[0];
            String chainIdentifier = chainId.split("_")[1].substring(0, 1);
            logger.info("handling {}",
                    chainId);
            Structure structure = StructureParser.source(pdbId).minimalParsing(true).parse();
            Optional<Chain> chain = structure.select()
                    .chainId(chainIdentifier)
                    .asOptionalChain();
            if (!chain.isPresent()) {
                logger.warn("{} is missing in structure",
                        chainId);
            }

            // check if OPM document can be retrieved
            Document opmDocument = OrientationsOfProteinsInMembranesAnnotator.getDocument(pdbId);

            // check for PLIP data
            Document plipDocument = PLIPIntraMolecularAnnotator.getDocument(chain.get());

            MembraneConstants.write(directory.resolve("pdb").resolve(pdbId + ".pdb"),
                    structure.getPdbRepresentation());
            MembraneConstants.write(directory.resolve("opm").resolve(pdbId + ".opm"),
                    opmDocument.html());
            MembraneConstants.write(directory.resolve("plip").resolve(chain.get().getChainIdentifier().getFullName() + ".plip"),
                    plipDocument.html());

            logger.info("{} is present",
                    chainId);
        } catch (ComputationException e) {
            logger.warn("no OPM or PLIP entry for {}",
                    chainId,
                    e);
        } catch (UncheckedIOException e) {
            logger.warn("biosciences.hs-mittweida.de seems to be unreachable");
        } catch (Exception e) {
            logger.warn("{} is obsolete",
                    chainId,
                    e);
        }
    }
}
