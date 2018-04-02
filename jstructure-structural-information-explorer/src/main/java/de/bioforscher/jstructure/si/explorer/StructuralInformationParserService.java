package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.explorer.model.ContactDistanceBin;
import de.bioforscher.jstructure.si.explorer.model.ContactStructuralInformation;
import de.bioforscher.jstructure.si.explorer.model.ResidueStructuralInformation;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class StructuralInformationParserService {
    private static final StructuralInformationParserService INSTANCE = new StructuralInformationParserService();

    private StructuralInformationParserService() {

    }

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
        Map<Pair<Integer, Integer>, List<String>> parsingMap = new HashMap<>();
        try(Stream<String> stream = new BufferedReader(new InputStreamReader(inputStream)).lines()) {
            stream.forEach(line -> {
                String[] split = line.split("\t");
                String[] idSplit = split[0].split(",");
                Pair<Integer, Integer> idPair = new Pair<>(Integer.valueOf(idSplit[0].split("\\(")[1].trim()),
                        Integer.valueOf(idSplit[1].split("\\)")[0].trim()));

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
                            .map(split -> Stream.of(split)
                                    .skip(8)
                                    .mapToDouble(Double::valueOf)
                                    .toArray())
                            .collect(Collectors.toList());

                    return new ContactStructuralInformation(entry.getKey().getLeft(),
                            entry.getKey().getRight(),
                            ContactDistanceBin.resolve(new Pair<>(IdentifierFactory.createResidueIdentifier(entry.getKey().getLeft()),
                                            IdentifierFactory.createResidueIdentifier(entry.getKey().getRight()))).orElse(null),
                            computeAverage(values, 0),
                            computeAverage(values, 1),
                            computeAverage(values, 2),
                            computeMaximum(values, 0),
                            computeMaximum(values, 1),
                            computeMaximum(values, 2),
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .map(ResidueIdentifier::getResidueNumber)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getLeft())) ||
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .map(ResidueIdentifier::getResidueNumber)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getRight())),
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .map(ResidueIdentifier::getResidueNumber)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getLeft())) &&
                            earlyFoldingResidues.stream()
                                    .map(Group::getResidueIdentifier)
                                    .map(ResidueIdentifier::getResidueNumber)
                                    .anyMatch(residueIdentifier -> residueIdentifier.equals(entry.getKey().getRight())));
                })
                // sort by average RMSD increase
//                .sorted(Comparator.comparingDouble(ContactStructuralInformation::getAverageRmsdIncrease).reversed())
                // sort by maximum RMSD increase
//                .sorted(Comparator.comparingDouble(ContactStructuralInformation::getMaximumRmsdIncrease).reversed())
                .collect(Collectors.toList());
    }

    private double computeAverage(List<double[]> values, int i) {
        return Double.valueOf(StandardFormat.format(values.stream()
                .mapToDouble(split -> split[i])
                .average()
                .orElse(0.0)));
    }

    private double computeMaximum(List<double[]> values, int i) {
        return Double.valueOf(StandardFormat.format(values.stream()
                .mapToDouble(split -> split[i])
                .max()
                .orElse(0.0)));
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
        int residueIdentifier = aminoAcid.getResidueIdentifier().getResidueNumber();
        List<ContactStructuralInformation> contactsOfAminoAcid = contacts.stream()
                .filter(contactStructuralInformation -> (contactStructuralInformation.getResidueIdentifier1() == residueIdentifier) ||
                        contactStructuralInformation.getResidueIdentifier2() == residueIdentifier)
                .collect(Collectors.toList());
        return new ResidueStructuralInformation(residueIdentifier,
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getAverageRmsdIncrease),
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getAverageTmScoreIncrease),
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getAverageQIncrease),
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getMaximumRmsdIncrease),
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getMaximumTmScoreIncrease),
                computeFeatureSum(contactsOfAminoAcid, ContactStructuralInformation::getMaximumQIncrease),
                earlyFoldingResidues.contains(aminoAcid));
    }

    private double computeFeatureSum(List<ContactStructuralInformation> contactsOfAminoAcid,
                                     ToDoubleFunction<ContactStructuralInformation> function) {
        return Double.valueOf(StandardFormat.format(contactsOfAminoAcid.stream()
                .mapToDouble(function)
                .sum()));
    }
}
