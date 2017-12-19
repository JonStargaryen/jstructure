package de.bioforscher.start2fold.collect;

import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Protein;
import org.codehaus.jackson.map.ObjectMapper;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Takes a snapshot of the Start2Fold database. Considers all chains to be named "A".
 */
public class A01_DownloadStart2FoldDataset {
    private static final Logger logger = LoggerFactory.getLogger(A01_DownloadStart2FoldDataset.class);

    public static void main(String[] args) throws IOException {
        Document proteinsDocument = Jsoup.connect(Start2FoldConstants.START2FOLD_PROTEINS_URL).get();

        List<Protein> proteins = proteinsDocument.getElementsByTag("tr").stream()
                // skip header
                .skip(1)
                .map(Protein::new)
                .collect(Collectors.toList());
        logger.info("current Start2Fold release contains {} proteins",
                proteins.size());

        // download XML entries
        for (Protein protein : proteins) {
            try {
                Start2FoldConstants.move(new URL(Start2FoldConstants.START2FOLD_URL + protein.getEntryId() + ".xml"),
                        Start2FoldConstants.XML_DIRECTORY.resolve(protein.getEntryId() + ".xml"));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        // download FASTA sequences
        // assumption: only first chain is considered - all chains are equal
        for (Protein protein : proteins) {
            try {
                String fastaSequence = Jsoup.connect(Start2FoldConstants.START2FOLD_URL + protein.getEntryId() + ".fasta").get()
                        .getElementsByTag("p")
                        .get(0)
                        .html()
                        .split("<br>")[1]
                        .trim();
                Start2FoldConstants.write(Start2FoldConstants.FASTA_DIRECTORY.resolve(protein.getEntryId() + ".fasta"),
                        ">" + protein.getEntryId() + System.lineSeparator() +
                        fastaSequence);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        // download PDB structures
        for (Protein protein : proteins) {
            try {
                Start2FoldConstants.move(new URL("https://files.rcsb.org/download/" + protein.getPdbId() + ".pdb"),
                        Start2FoldConstants.PDB_DIRECTORY.resolve(protein.getEntryId() + ".pdb"));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        // write JSON representation of dataset for easy mapping later on
        String jsonRepresentation = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(proteins);
        Start2FoldConstants.write(Start2FoldConstants.BASE_DIRECTORY.resolve("maintable.json"),
                jsonRepresentation);
    }
}
