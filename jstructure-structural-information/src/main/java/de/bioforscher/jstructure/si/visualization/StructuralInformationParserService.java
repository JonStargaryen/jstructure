package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class StructuralInformationParserService {
    private static final StructuralInformationParserService INSTANCE = new StructuralInformationParserService();

    public static StructuralInformationParserService getInstance() {
        return INSTANCE;
    }

    public List<ContactStructuralInformation> parseContactStructuralInformation(Path path,
                                                                                List<AminoAcid> earlyFoldingResidues) {
        try {
            return parseContactStructuralInformationFile(Files.newInputStream(path), earlyFoldingResidues);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public List<ContactStructuralInformation> parseContactStructuralInformationFile(InputStream inputStream,
                                                                                    List<AminoAcid> earlyFoldingResidues) {
        Map<Pair<ResidueIdentifier, ResidueIdentifier>, List<String>> parsingMap = new HashMap<>();
        try(Stream<String> stream = new BufferedReader(new InputStreamReader(inputStream)).lines()) {
            stream.forEach(line -> {
                String[] split = line.split("\t");
                String[] idSplit = split[0].split(",");
                ResidueIdentifier residueIdentifier1 = IdentifierFactory.createResidueIdentifier(idSplit[0].split("-")[1].trim());
                ResidueIdentifier residueIdentifier2 = IdentifierFactory.createResidueIdentifier(idSplit[1].split("-")[1].split("\\)")[0].trim());
                Pair<ResidueIdentifier, ResidueIdentifier> idPair = new Pair<>(residueIdentifier1, residueIdentifier2);

                if(!parsingMap.containsKey(idPair)) {
                    parsingMap.put(idPair, new ArrayList<>());
                }

                parsingMap.get(idPair).add(line);
            });
        }

        return parsingMap.entrySet()
                .stream()
                .map(entry -> {
                    List<double[]> values = entry.getValue()
                            .stream()
                            .map(line -> line.split("\t"))
                            .map(split -> new double[] {
                                    Double.valueOf(split[2]),
                                    Double.valueOf(split[3]),
                                    Double.valueOf(split[4])
                            })
                            .collect(Collectors.toList());

                    return new ContactStructuralInformation(entry.getKey().getLeft(),
                            entry.getKey().getRight(),
                            ContactDistanceBin.resolve(entry.getKey()).orElse(null),
                            values.stream()
                                    .mapToDouble(split -> split[0])
                                    .average()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> -split[1])
                                    .average()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> -split[2])
                                    .average()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> split[0])
                                    .max()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> -split[1])
                                    .max()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> -split[2])
                                    .max()
                                    .orElse(0.0),
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getLeft())) ||
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getRight())),
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getLeft())) &&
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getRight())));
                })
                // sort by average RMSD increase
                .sorted(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed())
                // sort by maximum RMSD increase
//                .sorted(Comparator.comparingDouble(ContactStructuralInformation::getMaximumRmsdIncrease).reversed())
                .collect(Collectors.toList());
    }

    public List<ResidueStructuralInformation> composeResidueStructuralInformation(List<AminoAcid> aminoAcids,
                                                                                  List<AminoAcid> earlyFoldingResidues,
                                                                                  List<ContactStructuralInformation> contacts) {
        return aminoAcids.stream()
                .map(aminoAcid -> composeResidueStructuralInformation(aminoAcid, earlyFoldingResidues, contacts))
                .collect(Collectors.toList());
    }

    private ResidueStructuralInformation composeResidueStructuralInformation(AminoAcid aminoAcid,
                                                                             List<AminoAcid> earlyFoldingResidues,
                                                                             List<ContactStructuralInformation> contacts) {
        ResidueIdentifier residueIdentifier = aminoAcid.getResidueIdentifier();
        List<ContactStructuralInformation> contactsOfAminoAcid = contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getResidueIdentifier1().equals(residueIdentifier) ||
                        contactStructuralInformation.getResidueIdentifier2().equals(residueIdentifier))
                .collect(Collectors.toList());
        return new ResidueStructuralInformation(residueIdentifier,
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getAverageRmsdIncrease)
                        .sum(),
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getAverageTmScoreIncrease)
                        .sum(),
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getAverageQIncrease)
                        .sum(),
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getMaximumRmsdIncrease)
                        .sum(),
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getMaximumTmScoreIncrease)
                        .sum(),
                contactsOfAminoAcid.stream()
                        .mapToDouble(ContactStructuralInformation::getMaximumQIncrease)
                        .sum(),
                earlyFoldingResidues.contains(aminoAcid));
    }
}
