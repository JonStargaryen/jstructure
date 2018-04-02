package de.bioforscher.jstructure.si.explorer.model;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.graph.ResidueGraph;
import de.bioforscher.jstructure.graph.ResidueGraphCalculations;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.explorer.StructuralInformationParserService;
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

    @Test
    public void shouldPrintStructuralInformationByContact() {
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

        ResidueGraph residueGraph = ResidueGraph.createDistanceResidueGraph(chain);
        ResidueGraphCalculations residueGraphCalculations = new ResidueGraphCalculations(residueGraph);
        contactStructuralInformation.stream()
                .map(si -> {
                    AminoAcid aminoAcid1 = chain.select()
                            .residueIdentifier(IdentifierFactory.createResidueIdentifier(si.getResidueIdentifier1()))
                            .asAminoAcid();
                    AminoAcid aminoAcid2 = chain.select()
                            .residueIdentifier(IdentifierFactory.createResidueIdentifier(si.getResidueIdentifier2()))
                            .asAminoAcid();
                    Pair<AminoAcid, AminoAcid> pair = new Pair<>(aminoAcid1, aminoAcid2);
                    double betweenness = residueGraphCalculations.betweenness(pair);
                    return StandardFormat.format(si.getAverageRmsdIncrease()) + "," +
                            StandardFormat.format(si.getAverageTmScoreIncrease()) + "," +
                            StandardFormat.format(si.getAverageQIncrease()) + "," +
                            StandardFormat.format(si.getMaximumRmsdIncrease()) + "," +
                            StandardFormat.format(si.getMaximumTmScoreIncrease()) + "," +
                            StandardFormat.format(si.getMaximumQIncrease()) + "," +
                            StandardFormat.format(betweenness);
                })
                .forEach(System.out::println);
    }

    @Test
    public void shouldPrintStructuralInformationByResidue() {
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
        List<ResidueStructuralInformation> residueStructuralInformation = StructuralInformationParserService.getInstance()
                .composeResidueStructuralInformation(aminoAcids,
                        earlyFoldingResidues,
                        contactStructuralInformation);

        ResidueGraph residueGraph = ResidueGraph.createDistanceResidueGraph(chain);
        ResidueGraphCalculations residueGraphCalculations = new ResidueGraphCalculations(residueGraph);
        residueStructuralInformation.stream()
                .map(si -> {
                    AminoAcid aminoAcid = chain.select()
                            .residueIdentifier(IdentifierFactory.createResidueIdentifier(si.getResidueIdentifier()))
                            .asAminoAcid();
                    double betweenness = residueGraphCalculations.betweenness(aminoAcid);
                    double closeness = residueGraphCalculations.closeness(aminoAcid);
                    double cc = residueGraphCalculations.clusteringCoefficient(aminoAcid);
                    int degree = residueGraph.degreeOf(aminoAcid);

                    return StandardFormat.format(si.getAverageRmsdIncrease()) + "," +
                            StandardFormat.format(si.getAverageTmScoreIncrease()) + "," +
                            StandardFormat.format(si.getAverageQIncrease()) + "," +
                            StandardFormat.format(si.getMaximumRmsdIncrease()) + "," +
                            StandardFormat.format(si.getMaximumTmScoreIncrease()) + "," +
                            StandardFormat.format(si.getMaximumQIncrease()) + "," +
                            StandardFormat.format(betweenness) + "," +
                            StandardFormat.format(closeness) + "," +
                            StandardFormat.format(cc) + "," +
                            StandardFormat.format(degree);
                })
                .forEach(System.out::println);
    }
}