package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class A06_EarlyFoldingFrequencies {
    private static int[] aminoAcidFrequencies = new int[20];
    private static int[] earlyFoldingFrequencies = new int[20];

    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.PANCSA_LIST)
                .forEach(A06_EarlyFoldingFrequencies::handleLine);

        System.out.println("aminoAcid,count,efr,count_freq,efr_freq");
        double total = IntStream.of(aminoAcidFrequencies)
                .sum();
        double efr = IntStream.of(earlyFoldingFrequencies)
                .sum();
        AminoAcid.Family.canonicalAminoAcids().forEach(aminoAcid -> {
            int ordinal = aminoAcid.ordinal();
            System.out.println(aminoAcid.getOneLetterCode() + "," +
                    StandardFormat.format(aminoAcidFrequencies[ordinal]) + "," +
                    StandardFormat.format(earlyFoldingFrequencies[ordinal]) + "," +
                    StandardFormat.format(aminoAcidFrequencies[ordinal] / total) + "," +
                    StandardFormat.format(earlyFoldingFrequencies[ordinal] / efr));
        });
    }

    private static void handleLine(String line) {
        System.out.println(line);
        String[] split = line.split(";");
        String entryId = split[0];
        String pdbId = split[1];
        List<Integer> experimentIds = Pattern.compile(",")
                .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                .map(Integer::valueOf)
                .collect(Collectors.toList());

        Structure structure = StructureParser.fromPdbId(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                experimentIds);

        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());

        chain.aminoAcids().forEach(aminoAcid -> {
            int ordinal = AminoAcid.Family.resolveOneLetterCode(aminoAcid.getOneLetterCode()).ordinal();
            aminoAcidFrequencies[ordinal]++;
        });

        earlyFoldingResidues.forEach(aminoAcid -> {
            int ordinal = AminoAcid.Family.resolveOneLetterCode(aminoAcid.getOneLetterCode()).ordinal();
            earlyFoldingFrequencies[ordinal]++;
        });
    }
}
