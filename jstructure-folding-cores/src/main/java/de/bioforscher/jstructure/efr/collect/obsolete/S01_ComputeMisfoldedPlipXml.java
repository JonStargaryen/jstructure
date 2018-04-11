package de.bioforscher.jstructure.efr.collect.obsolete;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.feature.interaction.PLIPRestServiceQuery;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class S01_ComputeMisfoldedPlipXml {
    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.DATA_DIRECTORY.resolve("native").resolve("wozniak.list"))
                .forEach(S01_ComputeMisfoldedPlipXml::handleMisfoldedProtein);

    }

    private static void handleMisfoldedProtein(String line) {
        try {
            String pdbId = line.split(";")[0].split("_")[0];
            System.out.println(pdbId);

            if(Files.exists(Paths.get("/home/bittrich/git/phd_sb_repo/data/native/xml/").resolve(pdbId + ".xml"))) {
                return;
            }

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.getFirstChain();

            Start2FoldConstants.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/native/xml/").resolve(pdbId + ".xml"),
                    PLIPRestServiceQuery.calculateIntraChainDocument(chain).html());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
