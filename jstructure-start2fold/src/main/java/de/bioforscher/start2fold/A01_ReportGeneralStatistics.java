package de.bioforscher.start2fold;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.model.FunctionalResidueAnnotation;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.FunctionalResidueParser;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import edu.northwestern.at.utils.math.statistics.FishersExactTest;
import org.jsoup.Jsoup;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Report general features of the dataset - checks for sanity.
 */
public class A01_ReportGeneralStatistics {
    private static List<Integer> early = new ArrayList<>();
    private static List<Integer> late = new ArrayList<>();
    private static List<Integer> functional = new ArrayList<>();
    private static List<Integer> nonFunctional = new ArrayList<>();
    private static List<Integer> overlap = new ArrayList<>();
    private static List<Integer> strong = new ArrayList<>();
    private static List<Integer> weak = new ArrayList<>();
    private static List<String> tableLines = new ArrayList<>();
    private static List<String> functionalTableLines = new ArrayList<>();
    private static int[] contingencyTable = new int[4];

    public static void main(String[] args) throws IOException {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());

        /*
         * Print general count statistics for early-folding-set defined by Pancsa et al., 2016.
         */
        Files.lines(Start2FoldConstants.PANCSA_LIST)
                .forEach(A01_ReportGeneralStatistics::handleEFRLine);
        int e = early.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int l = late.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int t = e + l;
        int f = functional.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int n = nonFunctional.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int ft = f + n;
        stringJoiner.add("EFR: " + early.size() + " proteins - " + t + " residues: " + e + " early, " + l + " late");
        stringJoiner.add("EFR-functional: " + functional.size() + " proteins - " + ft + " residues: " + f + " functional, " + n + " non-functional");
        stringJoiner.add("");

        /*
        Pancsa, 2016:
        For each of the 30 remaining proteins (Table S1 in the Sup-
        porting Material), totaling 3393 residues, the 482 resid ues that were
        experimentally determined as early folding were classified as belonging
        to class F, the remaining 2911 residues as belonging to class N. An
         */
        Files.list(Start2FoldConstants.XML_DIRECTORY)
                // ignore the file whose sequences cannot be aligned
                .filter(path -> !path.toFile().getName().contains("STF0034"))
                .forEach(A01_ReportGeneralStatistics::handleStrongFile);
        int s = strong.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int w = weak.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int t2 = s + w;
        stringJoiner.add("stability: " + strong.size() + " proteins - " + t2 + " residues: " + s + " strong, " + w + " weak");
        stringJoiner.add("");

        /*
         * Print contingency table between functional-early split.
         */
        stringJoiner.add("\\begin{tabular}{l r r }" + System.lineSeparator() +
                " & functional & non-functional \\\\ \\hline" + System.lineSeparator() +
                "early & " + contingencyTable[0] + " & " + contingencyTable[1] + " \\\\" + System.lineSeparator() +
                "late & " + contingencyTable[2] + " & " + contingencyTable[3] + "\\\\" + System.lineSeparator() +
                "\\end{tabular}");
        stringJoiner.add("");

        /*
         * Print functional test results.
         */
        stringJoiner.add(functionalTableLines.stream()
                .collect(Collectors.joining(System.lineSeparator(),
                        "\\begin{tabular}{ l r r r r r r }" + System.lineSeparator() +
                                "Start2Fold & residues & early & functional & overlap & \\textit{p}-value & significance \\\\ \\hline" + System.lineSeparator(),
                        " \\hline" + System.lineSeparator() +
                                "&&&& " + FishersExactTest.fishersExactTest(contingencyTable[0],
                                contingencyTable[1],
                                contingencyTable[2],
                                contingencyTable[3])[0] + " & ? \\\\" + System.lineSeparator() +
                                "\\end{tabular}")));
        stringJoiner.add("");

        /*
         * Print supplementary table for EFR.
         */
        tableLines.sort(Comparator.comparing(tableLine -> tableLine.split(" & ")[0]));
        stringJoiner.add(tableLines.stream()
                .collect(Collectors.joining(System.lineSeparator(),
                        "\\begin{tabular}{ l l l r r r r }" + System.lineSeparator() +
                        "Start2Fold & PDB & UniProt & residues & early & functional & overlap \\\\ \\hline" + System.lineSeparator(),
                        " \\hline" + System.lineSeparator() + "&&& " + t + " & " + e + " & " + f + " & " +
                                overlap.stream().mapToInt(Integer::valueOf).sum() + " \\\\" + System.lineSeparator() +
                                "\\end{tabular}")));

        System.out.println(stringJoiner.toString());
    }

    private static void handleStrongFile(Path path) {
        try {
            String pdbId = Jsoup.parse(path.toFile(), "UTF-8").getElementsByTag("protein").attr("pdb_id");
            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parse(chain, path);

            boolean hasStabilityData = chain.aminoAcids()
                    .map(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class))
                    .anyMatch(Start2FoldResidueAnnotation::isStrong);

            if(!hasStabilityData) {
                return;
            }

            long count = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isStrong())
                    .count();

            strong.add((int) count);
            weak.add((int) (chain.aminoAcids().count() - count));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void handleEFRLine(String line) {
        String[] split = line.split(";");
        String entryId = split[0];
        String pdbId = split[1];
        List<Integer> experimentIds = Pattern.compile(",")
                .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                .map(Integer::valueOf)
                .collect(Collectors.toList());
        int numberOfEarlyFoldingResidues = Integer.valueOf(split[3]);

        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                experimentIds);

        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());
        List<AminoAcid> lateFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> !earlyFoldingResidues.contains(aminoAcid))
                .collect(Collectors.toList());

        early.add(earlyFoldingResidues.size());
        late.add((int) (chain.aminoAcids().count() - earlyFoldingResidues.size()));

        if(earlyFoldingResidues.size() != numberOfEarlyFoldingResidues) {
            System.err.println("number of EFR did not match expectation for " + entryId + ": " +
                    earlyFoldingResidues.size() + " vs " + numberOfEarlyFoldingResidues);
        }

        String uniProtId = split[4];
        List<Integer> functionalResidueNumbers = Pattern.compile(",")
                .splitAsStream(split[5].replaceAll("\\[", "").replaceAll("]", ""))
                .flatMapToInt(numberString -> {
                    if(!numberString.contains("-")) {
                        return IntStream.of(Integer.valueOf(numberString));
                    }
                    String[] numberStringSplit = numberString.split("-");
                    return IntStream.rangeClosed(Integer.valueOf(numberStringSplit[0]),
                            Integer.valueOf(numberStringSplit[1]));
                })
                .boxed()
                .collect(Collectors.toList());
        List<AminoAcid> functionalResidues = new ArrayList<>();
        // do nothing if no annotation of functional residues exists
        if(!functionalResidueNumbers.isEmpty()) {
            FunctionalResidueParser.parse(chain, functionalResidueNumbers);
            chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(FunctionalResidueAnnotation.class).isFunctional())
                    .forEach(functionalResidues::add);
        }
        List<AminoAcid> nonFunctionalResidues = chain.aminoAcids()
                .filter(aminoAcid -> !functionalResidues.contains(aminoAcid))
                .collect(Collectors.toList());

        int earlyFunctionalCount = 0;
        if(!functionalResidues.isEmpty()) {
            functional.add(functionalResidues.size());
            nonFunctional.add((int) chain.aminoAcids().count() - functionalResidues.size());
            earlyFunctionalCount = SetOperations.intersection(earlyFoldingResidues, functionalResidues).size();
            overlap.add(earlyFunctionalCount);

            int ef = earlyFunctionalCount;
            int en = SetOperations.intersection(earlyFoldingResidues, nonFunctionalResidues).size();
            int lf = SetOperations.intersection(lateFoldingResidues, functionalResidues).size();
            int ln = SetOperations.intersection(lateFoldingResidues, nonFunctionalResidues).size();

            contingencyTable[0] += ef;
            contingencyTable[1] += en;
            contingencyTable[2] += lf;
            contingencyTable[3] += ln;

            double[] test = FishersExactTest.fishersExactTest(ef, en, lf, ln);
            System.out.println("values: " + ef + ", " + en + ", " + lf + ", " + ln);
            System.out.println("test: " + Arrays.toString(test));
            functionalTableLines.add(entryId + " & " +
                    chain.aminoAcids().count() + " & " +
                    earlyFoldingResidues.size() + " & " +
                    functionalResidues.size() + " & " +
                    ef + " & " +
                    StandardFormat.format(test[0]) + " & " +
                    "? \\\\");
        }

        tableLines.add(entryId + " & " +
                pdbId + "\\_A & " +
                uniProtId + " & " +
                chain.aminoAcids().count() + " & " +
                earlyFoldingResidues.size() + " & " +
                (functionalResidues.isEmpty() ? "-" : functionalResidues.size()) + " & " +
                (functionalResidues.isEmpty() ? "-" : earlyFunctionalCount) + " \\\\"
        );
    }
}
