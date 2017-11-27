package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

public class A01_DownloadPdbStructures {
    private static final Logger logger = LoggerFactory.getLogger(A01_DownloadPdbStructures.class);

    public static void main(String[] args) {
        List<String> pdbIds = ContactsConstants.list(ContactsConstants.START2FOLD_DIRECTORY.resolve("start2fold"))
                .flatMap(path -> {
                    try {
                        return Files.lines(path);
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                })
                .filter(line -> line.startsWith("#pdb: "))
                .map(line -> line.split("#pdb: ")[1])
                .map(line -> line.split("_")[0])
                .collect(Collectors.toList());

        ContactsConstants.write(ContactsConstants.START2FOLD_DIRECTORY.resolve("ids.list"),
                pdbIds.stream().collect(Collectors.joining(System.lineSeparator())));

        pdbIds.forEach(A01_DownloadPdbStructures::handlePdbId);
    }

    private static void handlePdbId(String pdbId) {
        try {
            logger.info("moving local pdb structure: {}", pdbId);
            Path pdbFile = Paths.get("/var/local/pdb/" + pdbId.substring(1, 3) + "/pdb" + pdbId + ".ent.gz");
            GZIPInputStream inputStream = new GZIPInputStream(Files.newInputStream(pdbFile));
            Path outDirectory = ContactsConstants.START2FOLD_DIRECTORY.resolve("pdb");
            Files.copy(inputStream, outDirectory.resolve(pdbId + ".pdb"));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
