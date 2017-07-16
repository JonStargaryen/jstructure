package studies.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPRestServiceQuery;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.nodes.Document;
import studies.StudyConstants;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Download OPM, PLIP and CIF definition files to speed up computations.
 * Created by bittrich on 6/7/17.
 */
public class ExternalResourceDownloader {
    public static void main(String[] args) throws IOException {
        MembraneConstants.PdbtmAlphaNr.getIds()
                .forEach(ExternalResourceDownloader::handleLine);
    }

    private static void handleLine(String line) {
        System.out.println(line);
        String pdbId = line.split("_")[0];

        // files may exist already
        if(MembraneConstants.PDBTM_OPM_PATH.resolve(pdbId + ".xml").toFile().exists()) {
            return;
        }

        try {
            Files.copy(new URL(String.format("https://files.rcsb.org/download/%s.pdb", pdbId)).openStream(), MembraneConstants.PDBTM_PDB_PATH.resolve(pdbId + ".pdb"));
        } catch (IOException e) {
            // thrown when multiple chains would result in the same file being created
            e.printStackTrace();
        }

        // load protein, calculate features
        Structure protein = StructureParser.source(pdbId).parse();

        // download ligand files
        protein.select()
                .ligands()
                .asFilteredGroups()
                .map(Group::getGroupPrototype)
                .map(GroupPrototype::getThreeLetterCode)
                .distinct()
                .forEach(ligandCode -> {
                    // do nothing when ligand files already exists
                    if (StudyConstants.list(MembraneConstants.PDBTM_CIF_PATH)
                            .map(Path::toFile)
                            .map(File::getName)
                            .anyMatch(name -> name.startsWith(ligandCode))) {
                        return;
                    }
                    Document cifDocument = GroupPrototypeParser.getDocument(ligandCode);
                    StudyConstants.write(MembraneConstants.PDBTM_CIF_PATH.resolve(ligandCode + ".xml"), cifDocument.html());
                });

        // fetch OPM files
        Document opmDocument = OrientationsOfProteinsInMembranesAnnotator.getDocument(protein.getProteinIdentifier().getPdbId());
        StudyConstants.write(MembraneConstants.PDBTM_OPM_PATH.resolve(protein.getProteinIdentifier().getPdbId() + ".xml"), opmDocument.html());

        // fetch PLIP results
        try {
            protein.chains().forEach(chain -> {
                Document plipDocument = PLIPRestServiceQuery.getDocument(chain.getChainIdentifier());
                StudyConstants.write(MembraneConstants.PDBTM_PLIP_PATH.resolve(chain.getChainIdentifier().getFullName() + ".xml"), plipDocument.html());
            });
        } catch (UncheckedIOException e) {
            e.getCause().printStackTrace();
        }
    }
}
