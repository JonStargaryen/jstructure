package de.bioforscher.jstructure.efr.analysis;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.EQuantScore;
import de.bioforscher.jstructure.efr.model.FunctionalResidueAnnotation;
import de.bioforscher.jstructure.efr.model.HotSpotScoring;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.EQuantParser;
import de.bioforscher.jstructure.efr.parser.EvolutionaryCouplingParser;
import de.bioforscher.jstructure.efr.parser.FunctionalResidueParser;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.geometry.GeometricProperties;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.graph.ResidueGraph;
import de.bioforscher.jstructure.graph.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A02_WriteDatasetCsv {
    private static final Logger logger = LoggerFactory.getLogger(A02_WriteDatasetCsv.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_NR_LIST)
                .map(A02_WriteDatasetCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "exposed,dCentroid," +
                                "terminusDistance," +
                                "plip_hydrogen,plip_hydrophobic,plip_bb,plip_total," +
                                "plip_l_hydrogen,plip_l_hydrophobic,plip_l_bb,plip_l_total," +
                                "plip_nl_hydrogen,plip_nl_hydrophobic,plip_nl_bb,plip_nl_total," +
                                "energy,egor,equant," +
                                "asa,loopFraction," +
                                "eccount,cumstrength,ecstrength,conservation," +
                                "plip_betweenness,plip_closeness,plip_clusteringcoefficient," +
                                "plip_hydrogen_betweenness,plip_hydrogen_closeness,plip_hydrogen_clusteringcoefficient," +
                                "plip_hydrophobic_betweenness,plip_hydrophobic_closeness,plip_hydrophobic_clusteringcoefficient," +
                                "conv_total,conv_l_total,conv_nl_total," +
                                "conv_betweenness,conv_closeness,conv_clusteringcoefficient," +
                                "plip_distinct_neighborhoods,conv_distinct_neighborhoods," +
                                "folds,functional,stable" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("foldingcores.csv"),
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

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parseStability(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"));
            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            EvolutionaryCouplingParser.parseHotSpotFile(chain,
                    Start2FoldConstants.COUPLING_DIRECTORY.resolve(entryId.toUpperCase() + "_hs.html"));

            EQuantParser.parseEQuantFile(chain,
                    Start2FoldConstants.EQUANT_DIRECTORY.resolve(entryId.toLowerCase() + ".equant-small.txt"));

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
            ResidueGraph conventionalProteinGraph = ResidueGraph.createResidueGraph(chain, ContactDefinitionFactory.createAlphaCarbonContactDefinition(8));

            return Optional.of(chain.aminoAcids()
                    .map(aminoAcid -> {
                        GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                        HotSpotScoring hotSpotScoring = aminoAcid.getFeature(HotSpotScoring.class);

                        PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                        PLIPInteractionContainer nonLocalPlipInteractionContainer = new PLIPInteractionContainer(null,
                                plipInteractionContainer
                                        .getInteractions()
                                        .stream()
                                        // interactions have to be non-local
                                        .filter(inter -> Math.abs(inter.getPartner1().getResidueIndex() - inter.getPartner2().getResidueIndex()) > 5)
                                        .collect(Collectors.toList()));
                        PLIPInteractionContainer localPlipInteractionContainer = new PLIPInteractionContainer(null,
                                plipInteractionContainer
                                        .getInteractions()
                                        .stream()
                                        // interactions have to be local
                                        .filter(inter -> !nonLocalPlipInteractionContainer.getInteractions().contains(inter))
                                        .collect(Collectors.toList()));

                        String equantScore = "NA";
                        try {
                            equantScore = StandardFormat.format(aminoAcid.getFeature(EQuantScore.class).getEvaluation());
                        } catch (ComputationException e) {
                            logger.warn("missing equant scoring for {}",
                                    aminoAcid);
                        }

                        String functionalAnnotation = "NA";
                        if(functionalResidues.size() > 0) {
                            functionalAnnotation = functionalResidues.contains(aminoAcid) ? "functional" : "non-functional";
                        }

                        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer = aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);
                        double terminusDistance = aminoAcids.indexOf(aminoAcid);
                        terminusDistance = Math.min(terminusDistance, aminoAcids.size() - terminusDistance);
                        terminusDistance /= (double) aminoAcids.size();

                        return pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +
                                sse.getSecondaryStructure().getReducedRepresentation() + "," +
                                sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                                sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                                (aminoAcid.getFeature(AccessibleSurfaceArea.class).isExposed() ? "exposed" : "buried") + "," +
                                StandardFormat.format(aminoAcid.getFeature(GeometricProperties.class).getDistanceToCentroid()) + "," +
                                StandardFormat.format(terminusDistance) + "," +

                                plipInteractionContainer.getHydrogenBonds().size() + "," +
                                plipInteractionContainer.getHydrophobicInteractions().size() + "," +
                                plipInteractionContainer.getBackboneInteractions().size() + "," +
                                plipInteractionContainer.getInteractions().size() + "," +

                                localPlipInteractionContainer.getHydrogenBonds().size() + "," +
                                localPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                                localPlipInteractionContainer.getBackboneInteractions().size() + "," +
                                localPlipInteractionContainer.getInteractions().size() + "," +

                                nonLocalPlipInteractionContainer.getHydrogenBonds().size() + "," +
                                nonLocalPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                                nonLocalPlipInteractionContainer.getBackboneInteractions().size() + "," +
                                nonLocalPlipInteractionContainer.getInteractions().size() + "," +

                                StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction()) + "," +
                                equantScore + "," +

                                StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," +

                                hotSpotScoring.getEcCount() + "," +
                                StandardFormat.format(hotSpotScoring.getCumStrength()) + "," +
                                StandardFormat.format(hotSpotScoring.getEcStrength()) + "," +
                                hotSpotScoring.getConservation() + "," +

                                StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getBetweenness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getCloseness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getClusteringCoefficient()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getBetweenness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getCloseness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrogenPlip().getClusteringCoefficient()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getBetweenness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getCloseness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getHydrophobicPlip().getClusteringCoefficient()) + "," +

                                conventionalProteinGraph.getContactsOf(aminoAcid).size() + "," +
                                conventionalProteinGraph.getLocalContactsOf(aminoAcid).size() + "," +
                                conventionalProteinGraph.getNonLocalContactsOf(aminoAcid).size() + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getBetweenness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getCloseness()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient()) + "," +

                                StandardFormat.format(residueTopologicPropertiesContainer.getFullPlip().getDistinctNeighborhoodCount()) + "," +
                                StandardFormat.format(residueTopologicPropertiesContainer.getConventional().getDistinctNeighborhoodCount()) + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                functionalAnnotation + "," +
                                (stableResidues.contains(aminoAcid) ? "stable" : "unstable");
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }
}
