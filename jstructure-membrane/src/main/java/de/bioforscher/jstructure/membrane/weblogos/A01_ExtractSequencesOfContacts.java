package de.bioforscher.jstructure.membrane.weblogos;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.stream.Collectors;

public class A01_ExtractSequencesOfContacts {
    public static void main(String[] args) {
        Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY
                .resolve("results");

        // tm contacts by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[32].equals("0"))
                .map(split -> new Pair<>(split[3], split[4]))
                .collect(Collectors.groupingBy(Pair::getLeft))
                .forEach((aminoAcid, value) -> {
                    String sequences = value
                            .stream()
                            .map(Pair::getRight)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-all")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[25].equals("0"))
                .map(split -> new Pair<>(split[3], split[4]))
                .collect(Collectors.groupingBy(Pair::getLeft))
                .forEach((aminoAcid, value) -> {
                    String sequences = value
                            .stream()
                            .map(Pair::getRight)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrogenbond")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[26].equals("0"))
                .map(split -> new Pair<>(split[3], split[4]))
                .collect(Collectors.groupingBy(Pair::getLeft))
                .forEach((aminoAcid, value) -> {
                    String sequences = value
                            .stream()
                            .map(Pair::getRight)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrophobic")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });

        // tm kinks by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have significant kink
                .filter(split -> split[11].equals("signficiant"))
                .map(split -> new Pair<>(split[3], split[4]))
                .collect(Collectors.groupingBy(Pair::getLeft))
                .forEach((aminoAcid, value) -> {
                    String sequences = value
                            .stream()
                            .map(Pair::getRight)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("kinks-all")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
    }
}
