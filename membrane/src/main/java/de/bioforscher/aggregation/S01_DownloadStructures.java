package de.bioforscher.aggregation;

import de.bioforscher.Constants;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Download all structures.
 * Created by bittrich on 2/13/17.
 */
public class S01_DownloadStructures {
    public static void main(String[] args) {
        Constants.lines(Paths.get(Constants.PDBTM_ALPHA_NR_LIST))
                .filter(Constants.isCommentLine.negate())
                // forget getChain identifier
                .map(line -> line.substring(0, 4))
                // ensure entries are unique
                .distinct()
                // download files and write them to directory
                .forEach(pdbId -> {
                    System.out.println(String.format("downloading '%s'", pdbId));
                    try (InputStream is = new URL(String.format(Constants.PDB_FETCH_URL, pdbId)).openStream()) {
                        Files.copy(is, Paths.get(Constants.STRUCTURE_PATH + pdbId + Constants.PDB_SUFFIX));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
    }
}
