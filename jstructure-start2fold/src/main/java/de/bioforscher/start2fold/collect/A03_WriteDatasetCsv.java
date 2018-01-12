package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.geometry.GeometricProperties;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.EQuantScore;
import de.bioforscher.start2fold.model.FunctionalResidueAnnotation;
import de.bioforscher.start2fold.model.HotSpotScoring;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.EQuantParser;
import de.bioforscher.start2fold.parser.EvolutionaryCouplingParser;
import de.bioforscher.start2fold.parser.FunctionalResidueParser;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class A03_WriteDatasetCsv {
    private static final Logger logger = LoggerFactory.getLogger(A03_WriteDatasetCsv.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A03_WriteDatasetCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "exposed,dCentroid," +
                                // PLIP contact definition
                                "plip_hydrogen,plip_hydrophobic,plip_bb,plip_total," +
                                "plip_l_hydrogen,plip_l_hydrophobic,plip_l_bb,plip_l_total," +
                                "plip_nl_hydrogen,plip_nl_hydrophobic,plip_nl_bb,plip_nl_total," +
                                "energy,egor,equant," +
                                "asa," +
                                "eccount,cumstrength,ecstrength,conservation," +
                                "folds,functional" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("foldingcores.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

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
                                        .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
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
                            equantScore = StandardFormat.format(aminoAcid.getFeature(EQuantScore.class).getEvaluation());                        } catch (ComputationException e) {
                            logger.warn("missing equant scoring for {}",
                                    aminoAcid);
                        }

                        String functionalAnnotation = "NA";
                        if(functionalResidues.size() > 0) {
                            functionalAnnotation = functionalResidues.contains(aminoAcid) ? "functional" : "non-functional";
                        }

                        return pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +
                                sse.getSecondaryStructure().getReducedRepresentation() + "," +
                                sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                                sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                                (aminoAcid.getFeature(AccessibleSurfaceArea.class).isExposed() ? "exposed" : "buried") + "," +
                                StandardFormat.format(aminoAcid.getFeature(GeometricProperties.class).getDistanceToCentroid()) + "," +

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

                                hotSpotScoring.getEcCount() + "," +
                                StandardFormat.format(hotSpotScoring.getCumStrength()) + "," +
                                StandardFormat.format(hotSpotScoring.getEcStrength()) + "," +
                                hotSpotScoring.getConservation() + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                                functionalAnnotation;
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
