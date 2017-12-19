package de.bioforscher.start2fold.collect;

import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Protein;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Either move old coupling files to new location or output FASTA sequences to compute missing scores.
 */
public class A02_RenameCouplingFiles {
    private static final Logger logger = LoggerFactory.getLogger(A02_RenameCouplingFiles.class);

    public static void main(String[] args) throws IOException {
        Document proteinsDocument = Jsoup.connect(Start2FoldConstants.START2FOLD_PROTEINS_URL).get();

        List<Protein> proteins = proteinsDocument.getElementsByTag("tr").stream()
                // skip header
                .skip(1)
                .map(Protein::new)
                .collect(Collectors.toList());
        logger.info("current Start2Fold release contains {} proteins",
                proteins.size());

        Path oldBasePath = Start2FoldConstants.DATA_DIRECTORY.resolve("datasets").resolve("foldingcores").resolve("couplings");

        for (Protein protein : proteins) {
            if(Files.exists(oldBasePath.resolve(protein.getPdbId() + "_A.fasta"))) {
                Start2FoldConstants.copy(oldBasePath.resolve(protein.getPdbId() + "_A.fasta"),
                        Start2FoldConstants.COUPLING_DIRECTORY.resolve(protein.getEntryId() + ".fasta"));
                Start2FoldConstants.copy(oldBasePath.resolve(protein.getPdbId() + "_A_ec.html"),
                        Start2FoldConstants.COUPLING_DIRECTORY.resolve(protein.getEntryId() + "_ec.html"));
                Start2FoldConstants.copy(oldBasePath.resolve(protein.getPdbId() + "_A_hs.html"),
                        Start2FoldConstants.COUPLING_DIRECTORY.resolve(protein.getEntryId() + "_hs.html"));
            } else {
                // move files to dedicated directory to compute coupling scores
                Start2FoldConstants.copy(Start2FoldConstants.FASTA_DIRECTORY.resolve(protein.getEntryId() + ".fasta"),
                        Paths.get("/home/bittrich/tmp/" + protein.getEntryId() + ".fasta"));
            }
        }
    }
}
