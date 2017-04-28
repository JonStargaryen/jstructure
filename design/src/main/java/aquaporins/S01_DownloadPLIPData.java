package aquaporins;

import de.bioforscher.jstructure.parser.plip.PLIPRestServiceQuery;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Fetch all PLIP data for the aquaporin data set.
 * Created by bittrich on 4/28/17.
 */
public class S01_DownloadPLIPData {
    public static void main(String[] args) {
        new SiftsMappingProvider().mapPfamIdToPdbIds("PF00230").forEach(pdbChainId -> {
            String document = PLIPRestServiceQuery.getPlipResults(pdbChainId.getPdbId().getPdbId(), pdbChainId.getChainId());
            try {
                Files.write(Paths.get(System.getProperty("user.home") + "/git/phd_sb_repo/data/aquaporins/plip/" + pdbChainId.getFullName() + ".xml"), document.getBytes());
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }
}
