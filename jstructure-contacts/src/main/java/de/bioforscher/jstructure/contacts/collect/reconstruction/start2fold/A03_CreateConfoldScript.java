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
//        ContactsConstants.write(DIRECTORY.resolve("plip-early-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForPlipEarly)
//                        .collect(Collectors.joining(System.lineSeparator())));
//        ContactsConstants.write(DIRECTORY.resolve("plip-sampled-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForPlipSampled)
//                        .collect(Collectors.joining(System.lineSeparator())));
//        ContactsConstants.write(DIRECTORY.resolve("plip-selected-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForPlipSelected)
//                        .collect(Collectors.joining(System.lineSeparator())));
//
//        ContactsConstants.write(DIRECTORY.resolve("conventional-early-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalEarly)
//                        .collect(Collectors.joining(System.lineSeparator())));
//        ContactsConstants.write(DIRECTORY.resolve("conventional-late-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalLate)
//                        .collect(Collectors.joining(System.lineSeparator())));
//        ContactsConstants.write(DIRECTORY.resolve("conventional-sampled-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalSampled)
//                        .collect(Collectors.joining(System.lineSeparator())));
//        ContactsConstants.write(DIRECTORY.resolve("conventional-selected-reconstruction.sh"),
//                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
//                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalSelected)
//                        .collect(Collectors.joining(System.lineSeparator())));
        ContactsConstants.write(DIRECTORY.resolve("conventional-late_sampled-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalLateSampled)
                        .collect(Collectors.joining(System.lineSeparator())));
        ContactsConstants.write(DIRECTORY.resolve("conventional-late_selected-reconstruction.sh"),
                ContactsConstants.lines(DIRECTORY.resolve("ids.list"))
                        .map(A03_CreateConfoldScript::handlePdbIdForConventionalLateSelected)
                        .collect(Collectors.joining(System.lineSeparator())));
    }

    private static String handlePdbIdForPlipEarly(String pdbId) {
        try {
            Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-early-plip-1/");
            Files.createDirectories(outDir);
            return "/home/bittrich/programs/confold_v1.0/confold.pl " +
                    "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                    "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-plip.rr " +
                    "-o " + outDir.toFile().getAbsolutePath();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForConventionalEarly(String pdbId) {
        try {
            Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-early-conventional-1/");
            Files.createDirectories(outDir);
            return "/home/bittrich/programs/confold_v1.0/confold.pl " +
                    "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                    "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-early-conventional.rr " +
                    "-o " + outDir.toFile().getAbsolutePath();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForConventionalLate(String pdbId) {
        try {
            Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-late-conventional-1/");
            Files.createDirectories(outDir);
            return "/home/bittrich/programs/confold_v1.0/confold.pl " +
                    "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                    "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-late-conventional.rr " +
                    "-o " + outDir.toFile().getAbsolutePath();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String handlePdbIdForPlipSampled(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-sampled-plip-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-plip-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForPlipSelected(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-selected-plip-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-selected-plip-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
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
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-sampled-conventional-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-conventional-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForConventionalLateSampled(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-late_sampled-conventional-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-late_sampled-conventional-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForConventionalSelected(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-selected-conventional-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-selected-conventional-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForConventionalLateSelected(String pdbId) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Path outDir = Paths.get("/home/bittrich/tmp/" + pdbId + "-late_selected-conventional-" + i + "/");
                Files.createDirectories(outDir);
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta " +
                        "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-late_selected-conventional-" + i + ".rr " +
                        "-o " + outDir.toFile().getAbsolutePath());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }
}
