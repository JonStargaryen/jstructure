package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class StructuralInformationParserServiceTest {
    @Test
    public void shouldParseStructuralInformationFile() {
        Chain chain = StructureParser.fromPdbId("1bdd")
                .parse()
                .getFirstChain();
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                TestUtils.getResourceAsInputStream("efr/STF0045.xml"),
                Stream.of(185).collect(Collectors.toList()));

        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());

        List<ContactStructuralInformation> contactStructuralInformation = StructuralInformationParserService.getInstance()
                .parseContactStructuralInformationFile(TestUtils.getResourceAsInputStream("si/STF0045.out"),
                        earlyFoldingResidues);

        contactStructuralInformation.stream()
                .map(ContactStructuralInformation::getCsvLine)
                .forEach(System.out::println);

        System.out.println();
        List<ResidueStructuralInformation> residueStructuralInformation = StructuralInformationParserService.getInstance()
                .composeResidueStructuralInformation(aminoAcids,
                        earlyFoldingResidues,
                        contactStructuralInformation);
        residueStructuralInformation.stream()
                .map(ResidueStructuralInformation::getCsvLine)
                .forEach(System.out::println);
    }
}