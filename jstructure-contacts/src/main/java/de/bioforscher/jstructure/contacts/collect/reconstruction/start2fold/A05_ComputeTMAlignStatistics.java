package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
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
                        "pdbId,mode,mapId,modelId,contacts,peraa,rmsd,tmscore" + System.lineSeparator(),
                        ""));

        ContactsConstants.write(ContactsConstants.START2FOLD_DIRECTORY.resolve("reconstruction-performance.csv"),
                output);
    }

    private static Optional<String> handleResultFile(Path tmalignFile) {
        String filename = tmalignFile.toFile().getName();
        String pdbId = filename.split("_")[0];
        String mode = tmalignFile.getParent().toFile().getName().split("-")[1];
        String mapId = "1";
        Path contactMapPath;
        if(mode.equals("naive")) {
            mapId = tmalignFile.getParent().toFile().getName().split("-")[2];
            contactMapPath = ContactsConstants.START2FOLD_DIRECTORY.resolve("maps")
                    .resolve(pdbId + "_A-naive-" + mapId + ".rr");
        } else {
            contactMapPath = ContactsConstants.START2FOLD_DIRECTORY.resolve("maps")
                    .resolve(pdbId + "_A-early-" + mode + ".rr");
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

        double rmsd = -1.0;
        double tmscore = -1.0;
        List<String> lines;
        try(Stream<String> stream = ContactsConstants.lines(tmalignFile)) {
            lines = stream.collect(Collectors.toList());
        }
        for(String line : lines) {
            if(line.startsWith("Aligned length=")) {
                rmsd = Double.valueOf(line.split("RMSD=")[1].trim().split(",")[0].trim());
            }
            if(line.startsWith("TM-score=")) {
                tmscore = Double.valueOf(line.split("TM-score=")[1].trim().split("\\(")[0].trim());
            }
        }

        if(rmsd == -1.0 || tmscore == -1.0) {
            return Optional.empty();
        }

        return Optional.of(pdbId + "," +
                mode + "," +
                mapId + "," +
                modelId + "," +
                numberOfContacts + "," +
                StandardFormat.format(numberOfContacts / (double) numberOfResidues) + "," +
                rmsd + "," +
                tmscore);
    }
}
