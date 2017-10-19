package de.bioforscher.jstructure.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
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
import java.util.stream.Collectors;

/**
 * Query the PDBTM and obtain the current list of non-redundant alpha-helical chain ids. Ensure that the corresponding
 * chain is present in the PDB, otherwise the chain is ignored and will not be present in the local list. Same is true
 * for the OPM and PLIP data. Everything ought to be sane for a chain to be added to the dataset.
 */
public class DatasetComposer {
    private static final Logger logger = LoggerFactory.getLogger(DatasetComposer.class);
    private final Path outputPath;

    public DatasetComposer(String fetchUrl, Path outputPath) {
        this.outputPath = outputPath;
        String selectedIds = MembraneConstants.lines(fetchUrl)
                .filter(this::handlePdbId)
                .collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(this.outputPath.resolve("ids.list"), selectedIds);
    }

    public DatasetComposer(Path outputPath) {
        this.outputPath = outputPath;
        String selectedIds = MembraneConstants.lines(outputPath.resolve("ids.list"))
                .filter(chainId -> handlePdbId(chainId, false))
                .collect(Collectors.joining(System.lineSeparator()));
        // overwrite list
        MembraneConstants.write(this.outputPath.resolve("ids.list"), selectedIds);
    }

    private boolean handlePdbId(String chainId, boolean parseOpm) {
        try {
            String pdbId = chainId.split("_")[0];
            String chainIdentifier = chainId.split("_")[1].substring(0, 1);
            logger.info("handling {} - parsing OPM: {}",
                    chainId,
                    parseOpm);
            Structure structure = StructureParser.source(pdbId).minimalParsing(true).parse();
            Optional<Chain> chain = structure.select()
                    .chainId(chainIdentifier)
                    .asOptionalChain();
            if (!chain.isPresent()) {
                logger.warn("{} is missing in structure",
                        chainId);
                return false;
            }

            // check if OPM document can be retrieved
            Document opmDocument = null;
            if(parseOpm) {
                OrientationsOfProteinsInMembranesAnnotator.getDocument(pdbId);
            }

            // check for PLIP data
            Document plipDocument = PLIPIntraMolecularAnnotator.getDocument(chain.get());

            MembraneConstants.write(outputPath.resolve("pdb").resolve(pdbId + ".pdb"),
                    structure.getPdbRepresentation());
            if(parseOpm) {
                MembraneConstants.write(outputPath.resolve("opm").resolve(pdbId + ".opm"),
                        opmDocument.html());
            }
            MembraneConstants.write(outputPath.resolve("plip").resolve(chain.get().getChainIdentifier().getFullName() + ".plip"),
                    plipDocument.html());

            logger.info("{} is present",
                    chainId);
            return true;
        } catch (ComputationException e) {
            logger.warn("no OPM or PLIP entry for {}",
                    chainId,
                    e);
        } catch (UncheckedIOException e) {
            logger.warn("biosciences.hs-mittweida.de seems to be unreachable");
            return false;
        } catch (Exception e) {
            logger.warn("{} is obsolete",
                    chainId,
                    e);
        }
        return false;
    }

    private boolean handlePdbId(String chainId) {
        return handlePdbId(chainId, false);
    }
}
