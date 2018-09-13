package de.bioforscher.jstructure.efr.analysis;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.FunctionalResidueAnnotation;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.FunctionalResidueParser;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A03_PrintStart2FoldDatasetTable {
    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.PANCSA_NR_LIST)
                .sorted(Comparator.comparing(line -> line.split(";")[0]))
                .map(A03_PrintStart2FoldDatasetTable::handleLine)
                .forEach(System.out::println);
    }

    private static String handleLine(String line) {
        try {
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parseStability(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"));
            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            List<AminoAcid> stableResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isStrong())
                    .collect(Collectors.toList());

            List<Integer> functionalResidueNumbers = Start2FoldConstants.extractFunctionalResidueNumbers(split);
            List<AminoAcid> functionalResidues = new ArrayList<>();
            // do nothing if no annotation of functional residues exists
            if(!functionalResidueNumbers.isEmpty()) {
                FunctionalResidueParser.parse(chain, functionalResidueNumbers);
                chain.aminoAcids()
                        .filter(aminoAcid -> aminoAcid.getFeature(FunctionalResidueAnnotation.class).isFunctional())
                        .forEach(functionalResidues::add);
            }

            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            long intersection = earlyFoldingResidues.stream()
                    .filter(functionalResidues::contains)
                    .count();

            return entryId + "\t" +
                    pdbId + "\t" +
                    split[2] + "\t" +
                    aminoAcids.size() + "\t" +
                    earlyFoldingResidues.size() + "\t" +
                    functionalResidues.size() + "\t" +
                    intersection;
        } catch (Exception e) {
            e.printStackTrace();
            return "";
        }
    }
}
