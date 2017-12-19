package de.bioforscher.start2fold;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.jsoup.Jsoup;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Report general features of the dataset - checks for sanity.
 */
public class A01_ReportGeneralStatistics {
    private static List<Integer> early = new ArrayList<>();
    private static List<Integer> late = new ArrayList<>();
    private static List<Integer> strong = new ArrayList<>();
    private static List<Integer> weak = new ArrayList<>();

    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.PANCSA_LIST)
                .forEach(A01_ReportGeneralStatistics::handleEFRLine);

        int e = early.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int l = late.stream()
                .mapToInt(Integer::valueOf)
                .sum();
        int t = e + l;
        System.out.println("EFR: " + early.size() + " proteins - " + t + " residues: " + e + " early, " + l + " late");

        /*
        Pancsa, 2016:
        For each of the 30 remaining proteins (Table S1 in the Sup-
        porting Material), totaling 3393 residues, the 482 residues that were
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
        System.out.println("stability: " + strong.size() + " proteins - " + t2 + " residues: " + s + " strong, " + w + " weak");
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

        long count = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .count();

        early.add((int) count);
        late.add((int) (chain.aminoAcids().count() - count));

        if(count != numberOfEarlyFoldingResidues) {
            System.err.println("number of EFR did not match expectation for " + entryId + ": " + count + " vs " + numberOfEarlyFoldingResidues);
        }
    }
}
