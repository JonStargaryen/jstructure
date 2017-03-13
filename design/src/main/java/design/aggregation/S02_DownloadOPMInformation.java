package design.aggregation;

import design.DesignConstants;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Fetches information on known pdb structures from the OPM.
 *
 * Diary:
 * <pre>
 *     Some (about the half) structures are missing as no *.opm information was available. Potentially, this is rather
 *     nice as we then have a data set to test resulting methods which is by no mean part of the training data set.
 *
 *     - 12/14/16 - 5A1S misses a ( which is crucial for parsing - add it when recreating the data
 *     - 03/13/17 - TODO fetch URL is flawed as some structures point to representative structure to avoid redundancy
 * </pre>
 * Created by S on 29.10.2016.
 */
@Deprecated
public class S02_DownloadOPMInformation {
    public static void main(String[] args) throws IOException {
        Files.list(Paths.get(DesignConstants.PDB_DIR))
             // extract id
             .map(path -> path.toFile().getName().substring(0, 4))
             // form OPM url
             .map(pdbId -> String.format(DesignConstants.OPM_FETCH_URL, pdbId))
             // fetch results
             .map(S02_DownloadOPMInformation::getDocument)
             // select missing entries
             .filter(S02_DownloadOPMInformation::filterAvailable)
             .forEach(document -> {
                 String pdbId = extractPdbId(document);
                 Elements data = document.select(".data");
                 try {
                     Files.write(Paths.get(DesignConstants.OPM_RAW_DIR + pdbId + ".opm"), data.toString().getBytes());
                 } catch (IOException e) {
                     e.printStackTrace();
                 }
             });
    }

    private static String extractPdbId(Document document) {
        return document.title().substring(0, 4);
    }

    private static boolean filterAvailable(Document document) {
        boolean present = !document.toString().contains("pdb/.pdb");
        if(!present) {
            System.err.println(extractPdbId(document) + " is missing");
        }
        return present;
    }

    private static Document getDocument(String url) {
        System.out.println("fetching " + url);
        try {
            return Jsoup.connect(url).get();
        } catch (IOException e) {
            e.printStackTrace();
            throw new UncheckedIOException(e);
        }
    }
}
