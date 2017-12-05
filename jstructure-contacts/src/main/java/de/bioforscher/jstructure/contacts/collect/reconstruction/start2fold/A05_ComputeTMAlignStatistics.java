package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A05_ComputeTMAlignStatistics {
    public static void main(String[] args) {
        String output = ContactsConstants.walk(ContactsConstants.START2FOLD_DIRECTORY.resolve("tmalign"))
                .filter(path -> !Files.isDirectory(path))
                .filter(path -> path.toFile().getName().endsWith(".tmalign"))
                .map(A05_ComputeTMAlignStatistics::handleResultFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdbId,strategy,scheme,mapId,modelId,contacts,peraa,rmsd,tmscore" + System.lineSeparator(),
                        ""));

        ContactsConstants.write(ContactsConstants.START2FOLD_DIRECTORY.resolve("reconstruction-performance.csv"),
                output);
    }

    private static Optional<String> handleResultFile(Path tmalignFile) {
        String filename = tmalignFile.toFile().getName();
        String pdbId = filename.split("_")[0];
        String strategy = tmalignFile.getParent().toFile().getName().split("-")[1];
        String scheme = tmalignFile.getParent().toFile().getName().split("-")[2];
        String mapId = "1";
        Path contactMapPath;
        System.out.println(pdbId + " > " + strategy + " > " + scheme + " > " + mapId);
        if(strategy.equals("early") || strategy.equals("full") || strategy.equals("deleted") || strategy.equals("late")) {
            contactMapPath = ContactsConstants.START2FOLD_DIRECTORY.resolve("maps")
                    .resolve(pdbId + "_A-" + strategy + "-" + scheme + ".rr");
        } else {
            mapId = tmalignFile.getParent().toFile().getName().split("-")[3];
            contactMapPath = ContactsConstants.START2FOLD_DIRECTORY.resolve("maps")
                    .resolve(pdbId + "_A-" + strategy + "-" + scheme + "-" + mapId + ".rr");
        }
        String modelId = tmalignFile.toFile().getName().split("model")[1].split("\\.")[0];
        int numberOfContacts;
        try(Stream<String> stream = ContactsConstants.lines(contactMapPath)) {
            numberOfContacts = (int) stream.count() - 1;
        }
        int numberOfResidues;
        try(Stream<String> stream = ContactsConstants.lines(ContactsConstants.START2FOLD_DIRECTORY.resolve("maps").resolve(pdbId + "_A.fasta"))) {
            numberOfResidues = stream.skip(1)
                    .mapToInt(String::length)
                    .findFirst()
                    .getAsInt();
        }

        try {
            ContactsConstants.TMAlignResult tmAlignResult = ContactsConstants.parseTMAlignResultFile(tmalignFile);
            return Optional.of(pdbId + "," +
                    strategy + "," +
                    scheme + "," +
                    mapId + "," +
                    modelId + "," +
                    numberOfContacts + "," +
                    StandardFormat.format(numberOfContacts / (double) numberOfResidues) + "," +
                    tmAlignResult.getRmsd() + "," +
                    tmAlignResult.getTmscore());
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
