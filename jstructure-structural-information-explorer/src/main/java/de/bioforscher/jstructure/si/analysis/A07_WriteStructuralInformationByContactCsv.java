package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.ContactDistanceBin;
import de.bioforscher.jstructure.efr.model.FunctionalResidueAnnotation;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.efr.parser.EvolutionaryCouplingParser;
import de.bioforscher.jstructure.efr.parser.FunctionalResidueParser;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.efr.parser.StructuralInformationParserService;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.interaction.PLIPInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.graph.ResidueGraph;
import de.bioforscher.jstructure.graph.ResidueGraphCalculations;
import de.bioforscher.jstructure.graph.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A07_WriteStructuralInformationByContactCsv {
    private static final Logger logger = LoggerFactory.getLogger(A07_WriteStructuralInformationByContactCsv.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("pancsa-si.list"))
                .filter(line -> line.startsWith("STF"))
                // skip entries lacking EFR annotation
                .filter(line -> !line.contains(";[];"))
                .map(A07_WriteStructuralInformationByContactCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res1,aa1,res2,aa2,distanceBin,long,medium,short,dCentroid," +
                                "avgRmsd,avgTm,avgQ,maxRmsd,maxTm,maxQ," +
                                "avgRmsdZ,numberOfTopScoringContacts," +
                                "plm,rank,topScoring02,topScoring04,topScoring06,topScoring08,topScoring10,topScoring12,topScoring14,topScoring16," +
                                "betweenness,avg_betweenness,avg_closeness,avg_clusteringcoefficient," +
                                "hydrogen,hydrophobic," +
                                "efr1,efr2,func1,func2,strong1,strong2,buried1,buried2,ordered1,ordered2,efrsse1,efrsse2,arom1,arom2," +
                                "efrAnnotation,strongAnnotation,functionalAnnotation,ecAnnotation" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("statistics").resolve("foldingcores-si-contacts.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

//            boolean sane = split[6].equalsIgnoreCase("true");

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra centroid = chain.calculate().centroid();

            Path start2foldXml = Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml");
            Start2FoldXmlParser.parseStability(chain,
                    start2foldXml);
            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    start2foldXml,
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

            List<AminoAcid> strongResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isStrong())
                    .collect(Collectors.toList());

            List<AminoAcid> orderedResidues = chain.aminoAcids()
                    .filter(aminoAcid -> !aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().isCoilType())
                    .collect(Collectors.toList());

            List<AminoAcid> buriedResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class).isBuried())
                    .collect(Collectors.toList());

            List<AminoAcid> residuesInEarlyFoldingSecondaryStructureElements = chain.aminoAcids()
                    .filter(aminoAcid -> !aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().isCoilType())
                    .filter(aminoAcid -> {
                        GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement = aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid);
                        List<AminoAcid> surroundingAminoAcids = chain.getAminoAcids().subList(surroundingSecondaryStructureElement.getStart(), surroundingSecondaryStructureElement.getEnd() + 1);
                        return surroundingAminoAcids.stream()
                                .anyMatch(earlyFoldingResidues::contains);
                    })
                    .collect(Collectors.toList());

            List<AminoAcid> aromaticResidues = chain.aminoAcids()
                    .filter(AminoAcid.Filter.AROMATIC)
                    .collect(Collectors.toList());

            List<ContactStructuralInformation> contactStructuralInformation = StructuralInformationParserService.getInstance()
                    .parseContactStructuralInformation(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("raw").resolve(entryId.toUpperCase() + ".out"),
                            chain,
                            earlyFoldingResidues);
            ResidueGraph conventionalProteinGraph = ResidueGraph.createResidueGraph(chain, ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0));
            ResidueGraphCalculations residueGraphCalculations = new ResidueGraphCalculations(conventionalProteinGraph);

            try {
                EvolutionaryCouplingParser.parsePlmScore(contactStructuralInformation,
                        Jsoup.parse(Start2FoldConstants.newInputStream(Start2FoldConstants.COUPLING_DIRECTORY.resolve(entryId + "_ec.html")), "UTF-8", ""),
                        chain.getAminoAcids().size());
            } catch (Exception e) {

            }
            boolean ecAnnotation = contactStructuralInformation.stream()
                    .anyMatch(csi -> csi.getPlmScore() != 0.0);

            PLIPInteractionContainer plipInteractionContainer = chain.getFeature(PLIPInteractionContainer.class);

            System.out.println("efr: " + (earlyFoldingResidues.size() > 0) + " strong: " + (strongResidues.size() > 0) + " functional: " + (functionalResidues.size() > 0) + " couplings: " + ecAnnotation);

            return Optional.of(contactStructuralInformation.stream()
                    .map(contact -> {
                        AminoAcid aminoAcid1 = chain.select()
                                .residueNumber(contact.getResidueIdentifier1())
                                .asAminoAcid();
                        AminoAcid aminoAcid2 = chain.select()
                                .residueNumber(contact.getResidueIdentifier2())
                                .asAminoAcid();
                        Pair<AminoAcid, AminoAcid> pair = new Pair<>(aminoAcid1, aminoAcid2);
                        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer1 = aminoAcid1.getFeature(ResidueTopologicPropertiesContainer.class);
                        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer2 = aminoAcid1.getFeature(ResidueTopologicPropertiesContainer.class);

                        LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra contactCentroid = aminoAcid1.calculate()
                                .centroid()
                                .add(aminoAcid2.calculate().centroid())
                                .divide(2);

                        return pdbId + "," +
                                "A" + "," +
                                contact.getResidueIdentifier1() + "," +
                                contact.getAa1() + "," +
                                contact.getResidueIdentifier2() + "," +
                                contact.getAa2() + "," +
                                contact.getContactDistanceBin() + "," +
                                (contact.getContactDistanceBin() == ContactDistanceBin.LONG) + "," +
                                (contact.getContactDistanceBin() == ContactDistanceBin.MEDIUM) + "," +
                                (contact.getContactDistanceBin() == ContactDistanceBin.SHORT) + "," +
                                StandardFormat.format(contactCentroid.distance(centroid)) + "," +

                                StandardFormat.format(contact.getAverageRmsdIncrease()) + "," +
                                StandardFormat.format(contact.getAverageTmScoreIncrease()) + "," +
                                StandardFormat.format(contact.getAverageQIncrease()) + "," +
                                StandardFormat.format(contact.getMaximumRmsdIncrease()) + "," +
                                StandardFormat.format(contact.getMaximumTmScoreIncrease()) + "," +
                                StandardFormat.format(contact.getMaximumQIncrease()) + "," +
                                StandardFormat.format(contact.getAverageRmsdIncreaseZScore()) + "," +
                                contact.getFractionOfTopScoringContacts() + "," +

                                StandardFormat.format(contact.getPlmScore()) + "," +
                                contact.getCouplingRank() + "," +
                                contact.istop02() + "," +
                                contact.isTop04() + "," +
                                contact.isTop06() + "," +
                                contact.isTop08() + "," +
                                contact.isTop10() + "," +
                                contact.isTop12() + "," +
                                contact.isTop14() + "," +
                                contact.isTop16() + "," +

                                StandardFormat.format(residueGraphCalculations.betweenness(pair)) + "," +
                                StandardFormat.format(0.5 * residueTopologicPropertiesContainer1.getConventional().getBetweenness() + 0.5 * residueTopologicPropertiesContainer2.getConventional().getBetweenness()) + "," +
                                StandardFormat.format(0.5 * residueTopologicPropertiesContainer1.getConventional().getCloseness() + 0.5 * residueTopologicPropertiesContainer2.getConventional().getCloseness()) + "," +
                                StandardFormat.format(0.5 * residueTopologicPropertiesContainer1.getConventional().getClusteringCoefficient() + 0.5 * residueTopologicPropertiesContainer2.getConventional().getClusteringCoefficient()) + "," +
                                //TODO potentially add average values of residue-level features

                                plipInteractionContainer.getHydrogenBonds()
                                        .stream()
                                        .anyMatch(hydrogenBond -> isContact(hydrogenBond, aminoAcid1, aminoAcid2)) + "," +
                                plipInteractionContainer.getHydrophobicInteractions()
                                        .stream()
                                        .anyMatch(hydrophobicInteraction -> isContact(hydrophobicInteraction, aminoAcid1, aminoAcid2)) + "," +

                                contact.isEarlyFoldingResidue() + "," +
                                contact.isEarlyFoldingContact() + "," +
                                residueIsInCollection(functionalResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(functionalResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                residueIsInCollection(strongResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(strongResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                residueIsInCollection(buriedResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(buriedResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                residueIsInCollection(orderedResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(orderedResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                residueIsInCollection(residuesInEarlyFoldingSecondaryStructureElements, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(residuesInEarlyFoldingSecondaryStructureElements, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                residueIsInCollection(aromaticResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                contactIsInCollection(aromaticResidues, contact.getResidueIdentifier1(), contact.getResidueIdentifier2()) + "," +
                                (earlyFoldingResidues.size() > 0) + "," +
                                (strongResidues.size() > 0) + "," +
                                (functionalResidues.size() > 0) + "," +
                                ecAnnotation;
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.info("calculation failed for {}\nby: {}",
                    line,
                    e.getMessage());
            return Optional.empty();
        }
    }

    private static boolean isContact(PLIPInteraction interaction, AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return (interaction.getPartner1().equals(aminoAcid1) && interaction.getPartner2().equals(aminoAcid2)) ||
                (interaction.getPartner2().equals(aminoAcid1) && interaction.getPartner1().equals(aminoAcid2));
    }

    private static boolean residueIsInCollection(List<AminoAcid> aminoAcids, int resNum1, int resNum2) {
        return aminoAcids.stream()
                .map(Group::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum1)) ||
                aminoAcids.stream()
                        .map(Group::getResidueIdentifier)
                        .map(ResidueIdentifier::getResidueNumber)
                        .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum2));
    }

    private static boolean contactIsInCollection(List<AminoAcid> aminoAcids, int resNum1, int resNum2) {
        return aminoAcids.stream()
                .map(Group::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum1)) &&
                aminoAcids.stream()
                        .map(Group::getResidueIdentifier)
                        .map(ResidueIdentifier::getResidueNumber)
                        .anyMatch(residueIdentifier -> residueIdentifier.equals(resNum2));
    }
}
