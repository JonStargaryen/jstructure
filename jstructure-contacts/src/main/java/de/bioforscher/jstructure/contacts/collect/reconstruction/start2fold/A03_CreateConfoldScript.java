package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A03_CreateConfoldScript {
    private static final Path DIRECTORY = ContactsConstants.START2FOLD_DIRECTORY;

    public static void main(String[] args) {
        String plipAverage = ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A03_CreateConfoldScript::handlePdbIdForAverage)
                .collect(Collectors.joining(System.lineSeparator()));
        ContactsConstants.write(DIRECTORY.resolve("plip-average-reconstruction.sh"), plipAverage);

        String plipMaximum = ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A03_CreateConfoldScript::handlePdbIdForMaximum)
                .collect(Collectors.joining(System.lineSeparator()));
        ContactsConstants.write(DIRECTORY.resolve("plip-maximum-reconstruction.sh"), plipMaximum);

        String sampledJobOutput = ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A03_CreateConfoldScript::handlePdbIdForSampled)
                .collect(Collectors.joining(System.lineSeparator()));
        ContactsConstants.write(DIRECTORY.resolve("sampled-reconstruction.sh"), sampledJobOutput);
    }

    private static String handlePdbIdForSampled(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-naive-" + i + "/"));
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-naive-" + i + ".rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-naive-" + i + "/");
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForAverage(String pdbId) {
        try {
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-avg/"));
            return "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-avg.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-avg/";
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForMaximum(String pdbId) {
        try {
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-max/"));
            return "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-max.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-max/";
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
