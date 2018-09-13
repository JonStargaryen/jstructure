package de.bioforscher.jstructure.efr.analysis;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.FunctionalResidueAnnotation;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.FunctionalResidueParser;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
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

public class A04_ComputeContactTypeFrequencies {
    private static int efr = 0;
    private static int func = 0;
    private static int efr_hb = 0;
    private static int efr_hi = 0;
    private static int func_hb = 0;
    private static int func_hi = 0;

    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.PANCSA_NR_LIST)
                .sorted(Comparator.comparing(line -> line.split(";")[0]))
                .forEach(A04_ComputeContactTypeFrequencies::handleLine);

        System.out.println(StandardFormat.format(100 * efr_hb / (double) efr));
        System.out.println(StandardFormat.format(100 * efr_hi / (double) efr));
        System.out.println(StandardFormat.format(100 * func_hb / (double) func));
        System.out.println(StandardFormat.format(100 * func_hi / (double) func));
    }

    private static void handleLine(String line) {
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

        List<Integer> functionalResidueNumbers = Start2FoldConstants.extractFunctionalResidueNumbers(split);
        List<AminoAcid> functionalResidues = new ArrayList<>();
        // do nothing if no annotation of functional residues exists
        if(!functionalResidueNumbers.isEmpty()) {
            FunctionalResidueParser.parse(chain, functionalResidueNumbers);
            chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(FunctionalResidueAnnotation.class).isFunctional())
                    .forEach(functionalResidues::add);
        }

        efr += earlyFoldingResidues.size();
        func += functionalResidues.size();

        efr_hb += earlyFoldingResidues.stream()
                .filter(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrogenBonds().size() > 0)
                .count();
        efr_hi += earlyFoldingResidues.stream()
                .filter(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrophobicInteractions().size() > 0)
                .count();

        func_hb += functionalResidues.stream()
                .filter(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrogenBonds().size() > 0)
                .count();
        func_hi += functionalResidues.stream()
                .filter(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrophobicInteractions().size() > 0)
                .count();
    }
}
