package design.aggregation;

import design.DesignConstants;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * Downloads pdb structures by a given list of ids.
 *
 * Diary:
 * <pre>
 *     Basis is the PDBTM nr alpha-helical id list.
 * </pre>
 * Created by S on 29.10.2016.
 */
public class S01_DownloadPDBStructures {
    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(DesignConstants.NR_ALPHA_IDS))
             .filter(DesignConstants.isCommentLine.negate())
             // forget chain identifier
             .map(line -> line.substring(0, 4))
             // ensure entries are unique
             .collect(Collectors.toSet())
             // download files and write them to directory
             .forEach(pdbId -> {
                 System.out.println(String.format("handling '%s'", pdbId));
                 try (InputStream is = new URL(String.format(DesignConstants.PDB_FETCH_URL, pdbId)).openStream()) {
                     Files.copy(is, Paths.get(DesignConstants.PDB_DIR + pdbId + DesignConstants.PDB_SUFFIX));
                 } catch (IOException e) {
                     e.printStackTrace();
                 }
             });
    }
}
