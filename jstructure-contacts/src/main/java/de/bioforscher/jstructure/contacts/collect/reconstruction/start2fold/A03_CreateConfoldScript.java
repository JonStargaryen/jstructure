package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringJoiner;
import java.util.stream.Collectors;

/**
 * FIXME: CONFOLD/CNS cannot handle long directory names - update script to write stuff to tmp directory and copy from there once finished
 */
@Deprecated
public class A03_CreateConfoldScript {
    private static final Path DIRECTORY = ContactsConstants.START2FOLD_DIRECTORY;

    public static void main(String[] args) {
        ContactsConstants.write(DIRECTORY.resolve("plip-early-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForPlipEarly)
                        .collect(Collectors.joining(System.lineSeparator())));
        ContactsConstants.write(DIRECTORY.resolve("plip-sampled-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForPlipSampled)
                        .collect(Collectors.joining(System.lineSeparator())));

        ContactsConstants.write(DIRECTORY.resolve("conventional-early-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalEarly)
                        .collect(Collectors.joining(System.lineSeparator())));
        ContactsConstants.write(DIRECTORY.resolve("conventional-sampled-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalSampled)
                        .collect(Collectors.joining(System.lineSeparator())));
    }

    private static String handlePdbIdForPlipEarly(String pdbId) {
        try {
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-early-plip-1/"));
            return "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-plip.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-early-plip-1/";
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForConventionalEarly(String pdbId) {
        try {
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-early-conventional-1/"));
            return "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-conventional.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-early-conventional-1/";
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForPlipSampled(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-sampled-plip-" + i + "/"));
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-plip-" + i + ".rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-sampled-plip-" + i + "/");
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForConventionalSampled(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-sampled-conventional-" + i + "/"));
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-conventional-" + i + ".rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/reconstructions/" + pdbId + "-sampled-conventional-" + i + "/");
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }
}
