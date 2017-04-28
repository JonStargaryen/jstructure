package aquaporins;

import de.bioforscher.jstructure.model.identifier.PdbChainId;
import de.bioforscher.jstructure.model.identifier.PdbId;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Downloads all PDB structures for the aquaporin data set.
 * Created by bittrich on 4/28/17.
 */
public class S02_DownloadStrctures {
    public static void main(String[] args) {
        new SiftsMappingProvider().mapPfamIdToPdbIds("PF00230").stream()
                .map(PdbChainId::getPdbId)
                .map(PdbId::getPdbId)
                .distinct()
                .forEach(pdbId -> {
                    try {
                        Files.copy(new URL("https://files.rcsb.org/download/" + pdbId + ".pdb").openStream(),
                            Paths.get(System.getProperty("user.home") + "/git/phd_sb_repo/data/aquaporins/structures/" + pdbId + ".pdb"));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
    }
}
