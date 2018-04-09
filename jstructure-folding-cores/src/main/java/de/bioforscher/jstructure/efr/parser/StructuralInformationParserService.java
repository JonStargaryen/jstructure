package de.bioforscher.jstructure.efr.parser;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.model.ContactDistanceBin;
import de.bioforscher.jstructure.efr.model.HotSpotScoring;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.efr.model.si.ReconstructionStructuralInformation;
import de.bioforscher.jstructure.efr.model.si.ResidueStructuralInformation;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

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

    public double[] parseAveragePerformance(InputStream inputStream) {
        try(Stream<String> stream = new BufferedReader(new InputStreamReader(inputStream)).lines()) {
            List<String> lines = stream.collect(Collectors.toList());
            double rmsd = lines.stream()
                    .map(line -> line.split("\t")[2])
                    .mapToDouble(Double::valueOf)
                    .average()
                    .orElse(0.0);
            double tm = lines.stream()
                    .map(line -> line.split("\t")[3])
                    .mapToDouble(Double::valueOf)
                    .average()
                    .orElse(0.0);
            double q = lines.stream()
                    .map(line -> line.split("\t")[4])
                    .mapToDouble(Double::valueOf)
                    .average()
                    .orElse(0.0);
            return new double[] { rmsd, tm, q };
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

        Chain chain = earlyFoldingResidues.get(0).getParentChain();
        Map<Pair<Integer, Integer>, List<ReconstructionStructuralInformation>> reconstructionMap = new HashMap<>();
        parsingMap.entrySet()
                .stream()
                .flatMap(entry -> {
                    String aa1 = chain.select()
                            .residueNumber(entry.getKey().getLeft())
                            .asAminoAcid()
                            .getOneLetterCode();
                    String aa2 = chain.select()
                            .residueNumber(entry.getKey().getRight())
                            .asAminoAcid()
                            .getOneLetterCode();

                    return entry.getValue()
                            .stream()
                            .map(line -> line.split("\t"))
                            .map(split -> new ReconstructionStructuralInformation(entry.getKey().getLeft(),
                                    aa1,
                                    entry.getKey().getRight(),
                                    aa2,
                                    ContactDistanceBin.resolve(new Pair<>(IdentifierFactory.createResidueIdentifier(entry.getKey().getLeft()),
                                            IdentifierFactory.createResidueIdentifier(entry.getKey().getRight()))).orElse(null),
                                    split[1].equals("true"),
                                    Double.valueOf(split[2]),
                                    Double.valueOf(split[3]),
                                    Double.valueOf(split[4]),
                                    Double.valueOf(split[5]),
                                    Double.valueOf(split[6]),
                                    Double.valueOf(split[7]),
                                    Double.valueOf(split[8]),
                                    Double.valueOf(split[9]),
                                    Double.valueOf(split[10])));
                })
                .forEach(rsi -> {
                    Pair<Integer, Integer> idPair = new Pair<>(rsi.getResidueIdentifier1(), rsi.getResidueIdentifier2());
                    if(!reconstructionMap.containsKey(idPair)) {
                        reconstructionMap.put(idPair, new ArrayList<>());
                    }

                    reconstructionMap.get(idPair).add(rsi);
                });

        List<ReconstructionStructuralInformation> reconstructionStructuralInformation = reconstructionMap.values()
                .stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
        int numberOfReconstructions = reconstructionStructuralInformation.size();
        double averageRmsd = reconstructionStructuralInformation.stream()
                .mapToDouble(ReconstructionStructuralInformation::getRmsdIncrease)
                .average()
                .orElse(0.0);
        double standardDeviationRmsd = new StandardDeviation().evaluate(reconstructionStructuralInformation.stream()
                .mapToDouble(ReconstructionStructuralInformation::getRmsdIncrease)
                .toArray());
        double averageMaximumRmsd = reconstructionMap.entrySet()
                .stream()
                .mapToDouble(entry -> entry.getValue().stream()
                        .mapToDouble(ReconstructionStructuralInformation::getRmsdIncrease)
                        .max()
                        .orElse(0.0))
                .average()
                .orElse(0.0);
        double standardDeviationMaximumRmsd = new StandardDeviation().evaluate(reconstructionMap.entrySet()
                .stream()
                .mapToDouble(entry -> entry.getValue().stream()
                        .mapToDouble(ReconstructionStructuralInformation::getRmsdIncrease)
                        .max()
                        .orElse(0.0))
                .toArray());

        List<ReconstructionStructuralInformation> topScoringReconstructions = reconstructionMap.values()
                .stream()
                .flatMap(Collection::stream)
                .sorted(Comparator.comparingDouble(ReconstructionStructuralInformation::getRmsdIncrease).reversed())
                .limit((int)  (0.1 * numberOfReconstructions))
                .collect(Collectors.toList());

        return reconstructionMap.entrySet()
                .stream()
                .map(entry -> {
                    List<ReconstructionStructuralInformation> values = entry.getValue();
                    ReconstructionStructuralInformation reference = values.get(0);

                    return new ContactStructuralInformation(reference.getResidueIdentifier1(),
                            reference.getAa1(),
                            reference.getResidueIdentifier2(),
                            reference.getAa2(),
                            reference.getContactDistanceBin(),
                            computeAverage(values, ReconstructionStructuralInformation::getBaselineRmsd),
                            computeAverage(values, ReconstructionStructuralInformation::getBaselineTmScore),
                            computeAverage(values, ReconstructionStructuralInformation::getBaselineQ),
                            computeAverage(values, ReconstructionStructuralInformation::getRmsdIncrease),
                            computeAverage(values, ReconstructionStructuralInformation::getTmScoreIncrease),
                            computeAverage(values, ReconstructionStructuralInformation::getqIncrease),
                            computeMaximum(values, ReconstructionStructuralInformation::getRmsdIncrease),
                            computeMaximum(values, ReconstructionStructuralInformation::getTmScoreIncrease),
                            computeMaximum(values, ReconstructionStructuralInformation::getqIncrease),
                            residueIsInCollection(earlyFoldingResidues, entry.getKey().getLeft(), entry.getKey().getRight()),
                            contactIsInCollection(earlyFoldingResidues, entry.getKey().getLeft(), entry.getKey().getRight()),
                            averageRmsd,
                            standardDeviationRmsd,
                            averageMaximumRmsd,
                            standardDeviationMaximumRmsd,
                            reconstructionStructuralInformation,
                            topScoringReconstructions);
                })
                .collect(Collectors.toList());
    }

    private boolean residueIsInCollection(List<AminoAcid> aminoAcids, int resNum1, int resNum2) {
        return aminoAcids.stream()
                .map(Group::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum1)) ||
                aminoAcids.stream()
                        .map(Group::getResidueIdentifier)
                        .map(ResidueIdentifier::getResidueNumber)
                        .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum2));
    }

    private boolean contactIsInCollection(List<AminoAcid> aminoAcids, int resNum1, int resNum2) {
        return aminoAcids.stream()
                .map(Group::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum1)) &&
                aminoAcids.stream()
                        .map(Group::getResidueIdentifier)
                        .map(ResidueIdentifier::getResidueNumber)
                        .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum2));
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

    private double computeAverage(List<ReconstructionStructuralInformation> reconstructionStructuralInformation,
                                  ToDoubleFunction<ReconstructionStructuralInformation> function) {
        return reconstructionStructuralInformation.stream()
                .mapToDouble(function)
                .average()
                .orElse(0.0);
    }

    private double computeMaximum(List<ReconstructionStructuralInformation> reconstructionStructuralInformation,
                                  ToDoubleFunction<ReconstructionStructuralInformation> function) {
        return reconstructionStructuralInformation.stream()
                .mapToDouble(function)
                .max()
                .orElse(0.0);
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
                aminoAcid.getOneLetterCode(),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getAverageRmsdIncrease),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getAverageTmScoreIncrease),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getAverageQIncrease),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getMaximumRmsdIncrease),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getMaximumTmScoreIncrease),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getMaximumQIncrease),
                earlyFoldingResidues.contains(aminoAcid),
                aminoAcid.getFeature(HotSpotScoring.class).getEcCount(),
                aminoAcid.getFeature(HotSpotScoring.class).getCumStrength(),
                aminoAcid.getFeature(HotSpotScoring.class).getConservation(),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getAverageRmsdIncreaseZScore),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getMaximumRmsdIncreaseZScore),
                computeFeatureAverage(contactsOfAminoAcid, ContactStructuralInformation::getFractionOfTopScoringContacts));
    }

    private double computeFeatureAverage(List<ContactStructuralInformation> contactOfAminoAcid,
                                         ToDoubleFunction<ContactStructuralInformation> function) {
        return Double.valueOf(StandardFormat.format(contactOfAminoAcid.stream()
                .mapToDouble(function)
                .average()
                .orElse(0.0)));
    }

    private double computeFeatureSum(List<ContactStructuralInformation> contactsOfAminoAcid,
                                     ToDoubleFunction<ContactStructuralInformation> function) {
        return Double.valueOf(StandardFormat.format(contactsOfAminoAcid.stream()
                .mapToDouble(function)
                .sum()));
    }
}
