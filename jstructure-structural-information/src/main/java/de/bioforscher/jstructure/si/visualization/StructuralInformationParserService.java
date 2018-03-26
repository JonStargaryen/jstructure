package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class StructuralInformationParserService {
    private static final StructuralInformationParserService INSTANCE = new StructuralInformationParserService();

    public static StructuralInformationParserService getInstance() {
        return INSTANCE;
    }

    public List<ContactStructuralInformation> parseStructuralInformationFile(Path path,
                                                                             List<AminoAcid> earlyFoldingResidues) {
        Map<Pair<ResidueIdentifier, ResidueIdentifier>, List<String>> parsingMap = new HashMap<>();
        try(Stream<String> stream = Files.lines(path)) {
            stream.forEach(line -> {
                String[] split = line.split("\t");
                String[] idSplit = split[0].split(",");
                ResidueIdentifier residueIdentifier1 = IdentifierFactory.createResidueIdentifier(idSplit[0].split("-")[1].trim());
                ResidueIdentifier residueIdentifier2 = IdentifierFactory.createResidueIdentifier(idSplit[1].split("-")[1].trim());
                Pair<ResidueIdentifier, ResidueIdentifier> idPair = new Pair<>(residueIdentifier1, residueIdentifier2);

                if(!parsingMap.containsKey(idPair)) {
                    parsingMap.put(idPair, new ArrayList<>());
                }

                parsingMap.get(idPair).add(line);
            });
        } catch (IOException e) {
            throw new UncheckedIOException(e);
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
                                    .mapToDouble(split -> split[1])
                                    .average()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> split[2])
                                    .average()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> split[0])
                                    .max()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> split[1])
                                    .max()
                                    .orElse(0.0),
                            values.stream()
                                    .mapToDouble(split -> split[2])
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
                .collect(Collectors.toList());
    }
}
