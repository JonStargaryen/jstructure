package de.bioforscher.jstructure.membrane.weblogos;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

        // top 10% tm contacts by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[32].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-all-t10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[25].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrogenbond-t10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[26].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrophobic-t10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });

        // bottom 10% tm contacts by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[32].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-all-b10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[25].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrogenbond-b10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[26].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom10(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrophobic-b10")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });

        // top 25% tm contacts by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[32].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-all-t25")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[25].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrogenbond-t25")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[26].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = top25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrophobic-t25")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });

        // bottom 25% tm contacts by aa
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[32].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-all-b25")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[25].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrogenbond-b25")
                                    .resolve(aminoAcid + ".txt"),
                            sequences);
                });
        MembraneConstants.lines(directory.resolve("aminoacids.csv"))
                .map(line -> line.split(";"))
                // transmembrane
                .filter(split -> split[8].equals("I"))
                // have helix-helix contact
                .filter(split -> !split[26].equals("0"))
                .map(AnnotatedAminoAcid::new)
                .collect(Collectors.groupingBy(AnnotatedAminoAcid::getAminoAcid))
                .forEach((aminoAcid, value) -> {
                    String sequences = bottom25(value)
                            .collect(Collectors.joining(System.lineSeparator()));
                    MembraneConstants.write(directory.resolve("weblogos")
                                    .resolve("contacts-hydrophobic-b25")
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

    private static Stream<String> bottom10(List<AnnotatedAminoAcid> annotatedAminoAcids) {
        int size = (int) (0.25 * annotatedAminoAcids.size());
        annotatedAminoAcids.sort(Comparator.comparingDouble(AnnotatedAminoAcid::getRigidity));
        return annotatedAminoAcids.stream()
                .limit(size)
                .map(AnnotatedAminoAcid::getSequence);
    }

    private static Stream<String> top10(List<AnnotatedAminoAcid> annotatedAminoAcids) {
        int size = (int) (0.25 * annotatedAminoAcids.size());
        annotatedAminoAcids.sort(Comparator.comparingDouble(AnnotatedAminoAcid::getRigidity).reversed());
        return annotatedAminoAcids.stream()
                .limit(size)
                .map(AnnotatedAminoAcid::getSequence);
    }

    private static Stream<String> bottom25(List<AnnotatedAminoAcid> annotatedAminoAcids) {
        int size = (int) (0.25 * annotatedAminoAcids.size());
        annotatedAminoAcids.sort(Comparator.comparingDouble(AnnotatedAminoAcid::getRigidity));
        return annotatedAminoAcids.stream()
                .limit(size)
                .map(AnnotatedAminoAcid::getSequence);
    }

    private static Stream<String> top25(List<AnnotatedAminoAcid> annotatedAminoAcids) {
        int size = (int) (0.25 * annotatedAminoAcids.size());
        annotatedAminoAcids.sort(Comparator.comparingDouble(AnnotatedAminoAcid::getRigidity).reversed());
        return annotatedAminoAcids.stream()
                .limit(size)
                .map(AnnotatedAminoAcid::getSequence);
    }

    static class AnnotatedAminoAcid {
        final String aminoAcid;
        final String sequence;
        final Double rigidity;

        AnnotatedAminoAcid(String[] split) {
            this.aminoAcid = split[3];
            this.sequence = split[4];
            this.rigidity = Double.valueOf(split[9]);
        }

        public String getAminoAcid() {
            return aminoAcid;
        }

        public String getSequence() {
            return sequence;
        }

        public Double getRigidity() {
            return rigidity;
        }
    }
}
